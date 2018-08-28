##################################################################
#                                                                #
#      MineClusters ALGORITHM: COMPUTATION OF CELL PATTERNS      #
#                                                                #
##################################################################

############################################
# Binary logistic regression               #
############################################

#' Perform binary logistic regression 
#'
#' Performs a binary logistic regression on the training data.
#' Also, if specified, performs validation tests.
#'
#' @param trainingSet list of level list x value list for training
#' @param validationSet list of level list x value list for validation
#' @return res list of model and accuracy index computed if validationSet is not NULL
#' 
#' @importFrom stats glm
#' @importFrom stats predict
#' @importFrom stats mean
#' @importFrom stats var
#' @importFrom stats binomial
#' 
#' @export
logistic_regression <- function(trainingSet, validationSet=NULL, print=F) {
	levels <- unique(trainingSet$levels)
	if (length(levels) == 1) levels <- c(levels, "FALSE")
	## To lessen the effects of multicollinearity                   ##
	## Scale "independent variables"                                ##
	trainValue <- (if (var(trainingSet$values)!=0) scale else identity)(trainingSet$values)
	trainLevel <- trainingSet$levels
	df <- data.frame(Level=trainLevel, Value=trainValue)
	model <- glm(Level ~ Value, data=df, family=binomial(link='logit'))
	if (!is.null(validationSet)) {
		validValue <- (if (var(validationSet$values)!=0) scale else identity)(validationSet$values)
		validLevel <- validationSet$levels
		nd <- data.frame(Level=validLevel, Value=validValue)
		results <- ifelse(predict(model, newdata=nd, type="response") > 0.5, levels[1], levels[2])
		res <- list(model=model, acc=1-mean(results != validLevel))
	}
	else res <- list(model=model, acc=NULL)
	if (!is.null(res$acc) && print) {
		clock(paste0("VALIDATION: Estimated accuracy for levels \'", levels[1], ", ", levels[2], 
			"\' is: ", round(res$acc, 3)))
	}
	return(res)
}

############################################
# Cross-validation: leave-one-out strategy #
############################################

#' Cross Validation function
#'
#' Performs cross validation on the data at step i.
#'
#' @param i integer identifier number of the training set
#' @param levels vector with categorical values
#' @param values vector with numeric values
#' @param sets set of possible training sets
#' @param f regression function with inputs training set and validation set
#' @return res list of model x associated accuracy coefficient obtained with set number i
#' 
#' @export
cross_validation_aux <- function(i, levels, values, sets, f, print=F) {
	p <- length(levels)
	validIDX <- (1:p)[!(1:p %in% sets[[i]])]
	trainLevels <- levels[sets[[i]]]
	trainValues <- values[sets[[i]]]
	validLevels <- levels[validIDX]
	validValues <- values[validIDX]
	res <- f(list(levels=trainLevels, values=trainValues), 
		list(levels=validLevels, values=validValues), print=print)
	return(res)
}

#' Perfom cross validation
#'
#' Performs cross validation on the data to determine the best model with regression.
#'
#' @param k integer: cutting the dataset in k parts
#' @param levels vector with categorical values
#' @param values vector with numeric values
#' @param f regression function with inputs training set and validation set
#' @param it current restart number
#' @param maxit maximum number of restart iterations
#' @param print logical for printing cross validation results
#' @return res list of model from the best performing training set x associated accuracy coefficient
#' 
#' @importFrom caret createFolds
#' 
#' @export
cross_validation <- function(k, levels, values, f=logistic_regression, it=0, maxit=1e1, print=F) {
	p <- length(levels)
	sets <- createFolds(1:p, k=k, list=T, returnTrain=F)
	bestMODEL <- NULL
	bestACC <- NULL
	tryCatch({
			fct <- function(i) return({cat(".", "", file=""); cross_validation_aux(i, levels, values, sets, f, print=print)})
			ls <- lapply(1:length(sets), fct)
			bestIDX <- which.max(sapply(ls, function(res) res$acc))
			bestACC <- ls[[bestIDX]]$acc
			bestMODEL <- ls[[bestIDX]]$model
		}, error=function(e) {
			clock(paste0("ERROR in pattern detection: {", e[[1]], "}\n"))
			if (it < maxit) {
				clock("Restart Cross Validation...\n")
				return(cross_validation(k, levels, values, it=it+1))
			}
			clock("Too many tries for Cross Validation... Restart procedure.\n")
			return(NULL)
		})
	if (!is.null(bestACC) && print) {
		clock(paste0("VALIDATION: Best (", k, "-fold cross-validation) accuracy for levels \'", 
			unique(levels)[1], ", ", unique(levels)[2], "\' is: ", round(bestACC, 3), "\n"))
	}
	return(list(model=bestMODEL, acc=bestACC))
}

############################################
# Compute expression patterns              #
############################################

#' Get pattern matrices and dropout/expressiveness logistic regression models 
#'
#' Gets pattern matrices and dropout/expressiveness logistic regression models .
#'
#' @param sca SingleCellAssay object that contains an aggregated (by replication)
#'            log-normalized expression matrix
#' @return res list with dropout and expressiveness models 
#'		and pattern matrix and pattern-oriented threshold matrix
#' 
#' @importFrom SummarizedExperiment assay
#' @importFrom MAST thresholdSCRNACountMatrix
#' 
#' @export
getPatterns_aux <- function(sca, nbins=floor(nrow(sca)/1e2), min_per_bin=2, fast=F) {
	p <- nrow(assay(sca))
	m <- ncol(assay(sca))
	capture.output(res <- thresholdSCRNACountMatrix(assay(sca), nbins=nbins, min_per_bin=min_per_bin))
	## Dropout values put at zero (baseline) value                                           ## 
	assay(sca) <- res$counts_threshold
  	## Get the gene expression thresholds                                                    ##
	threshold <- sapply(as.character(res$bin), function(bin) res$cutpoint[bin])
	values_mild <- assay(sca)
	levels_mild <- apply(values_mild, 2, function(x) ifelse(x > threshold, "high", "mild"))
	## Resulting pattern matrix                                                              ##
	## P_ij = kronecker(M_ij > threshold[j])                                                 ##
	P <- apply(levels_mild == "high", 2, as.integer)
	if (fast) return(list(dropout=NULL, mild=NULL, pattern=P, threshold=NULL, peaks=NULL))
	## Resulting pattern-oriented threshold matrix                                           ##
	Q <- matrix(sapply(1:m, function(j) ifelse(P[, j] == 1, -1, 1)*threshold[j]), byrow=F, nrow=p)
	colnames(Q) <- colnames(P)
	clock("(2)a: Pattern matrices")
	## Dropout events depend on the considered cell                                          ##
	## and the expression value of the gene across the cells                                 ##
	values_dropout <- rowSums(apply(assay(sca), 2, function(x) if (sum(x) == 0) rep(0, p) else x/sum(x)))
	levels_dropout <- ifelse(values_dropout == 0, "dropout", "not_dropout")
  	## Matrix with column "Dropout" and column "Median" (log2(TPM+1) value across all cells) ##
	## Performs binary logistic regression to fit dropout distribution                       ##
	model_dropout <- cross_validation(k=8, levels_dropout, values_dropout)$model
	if (is.null(model_dropout)) return(NULL)
	clock("(2)b: Dropout model fitted")
	## Matrix with column "Mild" and column "Value" (log2(TPM+1) value in each cell)         ##
	## Performs binary logistic regression to fit expression distribution for each condition ##
	model_mild <- lapply(1:m, function(j){cat("*", "", file="");cross_validation(k=5, levels_mild[, j], unlist(values_mild[, j]))$model})
	if (any(sapply(model_mild, is.null))) return(NULL)
	clock("(2)c: Expressiveness model fitted")
	## Get the two highest peaks, if they exist, for each cell                               ##
	peaks <- lapply(as.character(res$bin), function(bin) {
		if (length(res$peaks[bin][[1]]) == 1) rep(res$peaks[bin][[1]][1, 1], 2)
		else {
			if (length(res$peaks[bin][[1]]) == 0) rep(NA, 2) 
			else res$peaks[bin][[1]][1:2, 1]
			}
	})
	return(list(dropout=model_dropout, mild=model_mild, pattern=P, threshold=Q, peaks=peaks))
}

getPatterns <- function(sca, options, nbinsden=1e2) {
	msgThreshold="No bimodal bins.  Try decreasing `min_per_bin` and/or increasing `num_bins`."
	tryCatch(
		return(getPatterns_aux(sca, nbins=nrow(sca)/nbinsden, fast=options$fast)), 
		error=function(e) {
			if (nbinsden > 1 && e[[1]] == msgThreshold) return(getPatterns(sca, nbinsden/10, fast=options$fast))
			else {
				clock(paste0("ERROR: ", e[[1]], "\n"))
				return(NULL)
			}
		})
}
