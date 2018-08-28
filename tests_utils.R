####################################################################
#                                                                  #
#               TEST file                                          #
#                                                                  #
####################################################################

## TOOLS            ##

verboseIt <- function(msg) print(paste0("MSG: ", msg, "!"))

formatParam <- function(param) {
	paste0(sapply(1:length(param), function(i) {
		paste0(sapply(1:length(param[i]), function(j) {
			paste0(names(param)[i], " = ", 
				paste0(as.vector(param[i][[j]]), collapse="|"))
			}),
		collapse=" ") 
	}), collapse="; ")
}

generateInterval <- function(step, start, end, nmodes=2) {
	return(list(interval=seq(start, end, step), 
		intervalprint=paste0(rep(paste0("[", start, ",", end, "]"), 
			nmodes), collapse="x")))
}

generateVectors <- function(l, nmodes=2) {
	return(lapply(apply(do.call(expand.grid, 
		rep(list(l$interval), nmodes)), 1, list), unlist))
}

heatmap_c <- function(mat, ordered=F, title="Heatmap") {
	require(gplots)
	if (ordered) {
		r <- order(sapply(colnames(mat), function(e) strsplit(e, "_")[[1]][2]))
		return(heatmap.2(mat[r, r], col=RColorBrewer::brewer.pal(9,"Blues"),
			revC=T, symm=T, Rowv=NULL, trace="none",
			main=title))
	}
	return(heatmap.2(mat, col=RColorBrewer::brewer.pal(9,"Blues"),
		revC=T, symm=T, trace="none",
		main=title))
}

######################
#    UTILS           #
######################

source("import.R")

#____________________#
# Computation utils  #
#____________________#

## f is a continuous nondecreasing function in [xmin, xmax]
bisectionInv <- function(f, x, eps=1e-3, step=10, lim=1e-2, maxiter=10, xmin=0, xmax=1, limInf=1e-2, default=NULL) {
	## Special cases                          ##
	if (x > xmax | x < xmin) return(-1)
	if (x == xmax) return(Inf)
	if (x == xmin) return(-Inf)
	## Generate search space                  ##
	l <- 0
	r <- step
	while((x-f(r) > eps) && maxiter > 0) {r <- r+step; maxiter <- maxiter-1}
	if (maxiter == 0) {
		#print(paste0("x=", round(x, 3), " r=", r, " diff=", round(x-f(r), 3)))
		#plot(seq(0, r, 1e-2), sapply(seq(0, r, 1e-2), f), type="l", col="blue")
		## If it fails, then returns "random" value ##
		return(default)
	}
	## While the search space is large enough ##
	while(r-l > lim) {
		m <- (r+l)/2
		## If f(m) == x %numeric approximation in R ##
		if (abs(f(m)-x) < eps) {
			## Inferior value ##
			while(abs(f(m)-x) < eps) m <- m-limInf
			return(m+limInf)
		}
		if (f(m) > x) r <- m else l <- m
	}
	return(m)
}

pmargin <- function(q, params, i, log.p=FALSE, lower.tail=TRUE) {
	## Same margin parameters for all copulas ##
	margin <- params[[1]]$paramsMargins[[i]]
	y <- margin$alpha*pgamma(q, shape=margin$lambda, scale=1)
		+(1-margin$alpha)*pgamma(q, shape=margin$mu/(1+margin$delta*margin$mu), 
			scale=1+margin$mu*margin$delta)
	if (lower.tail) return(ifelse(log.p, log(y), y))
	else return(ifelse(log.p, log(1-y), 1-y))
}

qmargin <- function(p, params, lower.tail=TRUE, log.p=FALSE) {
	## Same margin parameters for all copulas ##
	margin <- params[[1]]$paramsMargins[[i]]
	return(bisectionInv(function(q) return(sapply(q, function(x) pmargin(ifelse(lower.tail, x, 1-x), params, i))), 
		ifelse(log.p, exp(p), p), default=sample(c(mu, 0), 1)/sample(1:2, 1)))
}

transformQ <- function(x, params) {
	q <- qnorm(as.vector(apply(as.matrix(1:length(x), ncol=1), 1, function(i) {
		pmargin(abs(x[i]), params, i, lower.tail=!(x[i]<0))
	})))
	return(q)
}

copula <- function(q, params, j) {
	param <- params[[j]]
	retvals <- sum(sapply(1:length(param$mixCop), function(k) {
		param$mixCop[[k]]*pmnormCustom(x=as.numeric(q), mean=param$muCop[[k]], 
		varcov=param$rhoCop[[k]], utrho=param$utrhoCop[[k]])
	}))
	return(retvals)
}

## e.g. params <- list(muCop=rep(list(1:2), 2), mixCop=list(0.5, 0.5), 
## 		rhoCop=list(diag(1, 2), diag(1, 2)), 
## 		paramsMargins=list(list(alpha=0.3, lambda=0.1, mu=4, delta=0.4), 
## 				list(alpha=0.4, lambda=0.3, mu=3, delta=0.1)))
rCopula <- function(params, n) {
	require(mnormt)
	nmodes <- length(params$mixCop)
	m <- sample(1:nmodes, prob=unlist(params$mixCop), size=n, replace=T)
	rvectors <- t(matrix(unlist(lapply(m, function(k) rmnorm(n=1, 
		mean=params$muCop[[k]], 
		varcov=params$rhoCop[[k]]))), nrow=n))
	pvectors <- pnorm(rvectors)
	p <- nrow(pvectors)
	uvectors <- do.call(cbind, lapply(1:n, function(j) sapply(1:p, function(i) {
		qmargin(pvectors[i, j], list(params))
	})))
	return(as.matrix(uvectors, nrow=d, ncol=n))
}

upperTriCustom <- function(x) return(.Internal(row(dim(x))) < .Internal(col(dim(x))))

#' Tests copula distribution inference
#'
#' @param d dimension of copulas
#' @param nmodes number of bimodal copulas
#' @param size size of sample matrix
#' 
#' FUNCTION: `mineClusters`
#' FILE: `mineClusters.R`
test_computation <- function(d=10, nmodes=2, size=5) {
	source("mineClusters.R")
	require(mclust)
	require(Matrix)
	verboseIt("-- START TEST")
	sp <- function() sample(seq(0, 1, 0.1), 1)
	## Same margins            ##
	paramsMargins <- lapply(1:d, function(i) list(alpha=sp(), lambda=0.1, mu=sp(), delta=sp()))
	## Same correlation matrix ##
	rho <- nearPD(matrix(rgamma(d*d, shape=1e-3, scale=1), d, d), corr=T)$mat
	rhoCop <- rep(list(rho), 2)
	utrhoCop <- rep(list(upperTriCustom(rho)), 2)
	params <- lapply(1:nmodes, function(mode) {
		muCop <- list(rep(0, d)+mode*1e-1, rep(log(mode+1), d))
		mix <- sp()
		mixCop <- list(mix, 1-mix)
		list(muCop=muCop, mixCop=mixCop, rhoCop=rhoCop, utrhoCop=utrhoCop, paramsMargins=paramsMargins)
	})
	## Creation of a matrix   ##
	M <- matrix(unlist(lapply(1:nmodes, function(i) {
			print(paste0("Create mode #", i))
			rCopula(params[[i]], size)
	})), ncol=nmodes, nrow=d, byrow=F)
	## "True" labels          ##
	vcell <- as.vector(as.matrix(sapply(1:nmodes, function(i) rep(i, size)), nrow=1))
	colnames(M) <- sapply(1:ncol(M), function(e) paste0("s#", e, "_cell", vcell[e]))
	rownames(M) <- sapply(1:nrow(M), function(e) paste0("feature_", e))
	names(vcell) <- sapply(1:ncol(M), function(e) paste0("s#", e))
	if (any(M<0)) {print("Matrix is not positive"); return(NULL)}
	## Launch mineClusters    ##
	tryCatch({
		res <- mineClusters(M, list(replicates=NULL, X=1, thres_fs=0.9, fast=F, parallel=T, normalized=F, islog=F))
	}, error=print)
	print("RESTART FAST VERSION")
	capture.output(resfast <- mineClusters(M, list(replicates=NULL, X=1, thres_fs=0.9, fast=T, parallel=F, normalized=F, islog=F)))
	## Results                ##
	#heatmap_c(res$dissimilarity, ordered=T)
	if (!is.null(res$clusters)) {
		print(paste0("ARI for regular version = ", round(adjustedRandIndex(res$clusters, vcell), 3)))
	}
	print(paste0("ARI for fast version = ", round(adjustedRandIndex(resfast$clusters, vcell), 3)))
	if (nmodes==2) {
		clusters <- normalize_labels(if (is.null(res$clusters)) resfast$clusters else res$clusters)
		#names(clusters) <- names(vcell)
		#cat("Clustering\n", "", file="")
		#print(clusters)
		#cat("True labels\n", "", file="")
		#print(vcell)
		cp <- function(a, b=a) return(sum(sapply(1:length(vcell), function(i) clusters[i]==a&&vcell[i]==b)))
		confusion <- matrix(c(cp(1), cp(2, 1), cp(1, 2), cp(2)), ncol=2, nrow=2, byrow=T)
		colnames(confusion) <- c("Predicted Class 1", "Predicted Class 2")
		rownames(confusion) <- c("True Class 1", "True Class 2")
		print(confusion)
	}
	verboseIt("-- END TEST")
	return(list(init=list(params=params), copulas=res$copulas, dissimilarity=res$dissimilarity, mat=res$matrix))
}

#' Plots true and infered copulas
#' @param res output of function `test_computation`
#' @param i integer identifier of the copula that should be plotted
plotCopulas <- function(res, i) {
	ee <- seq(0, max(res$mat), 0.1)
	x <- lapply(ee, function(e) rep(e, nrow(res$mat)))
	print("Compute true values...")
	y <- sapply(x, function(xi) copula(transformQ(xi, res$init$params), res$init$params, i)
	print("Compute infered values...")
	yy <- sapply(x, function(xi) copula(transformQ(xi, res$copulas$params), res$copulas$params, i)
	plot(ee, y, type="l", col="blue", title="Comparison of copulas: true in blue, infered in red", ylim=c(-0.1, 1),
		xlab="Values of coordinates", ylab="Cumulative Probability")
	lines(ee, yy, type="l", col="red")
	return(list(x=x, initial=y, infered=y))
}

#____________________#
# Pattern utils      #
#____________________#

#' Tests pattern inference part
#'
#' @param model a R character string for the model located in ../datasets/
#' 
#' FUNCTION: `logistic_regression`
#' FILE: `pattern.R`
#'
#' FUNCTION: `cross_validation`
#' FILE: `pattern.R`
#'
#' FUNCTION: `filterMatrix`
#' FILE: `feature_selection.R`
#'
#' FUNCTION: `normalizeLogData`
#' FILE: `feature_selection.R`
#'
#' FUNCTION: `matrix2sca`
#' FILE: `feature_selection.R`
#'
#' FUNCTION: `filterGE`
#' FILE: `feature_selection.R`
#'
#' FUNCTION: `featureSelection`
#' FILE: `feature_selection.R`
#'
#' FUNCTION: `group_by_condition`
#' FILE: `utils.R`
#'
#' FUNCTION: `getPatterns`
#' FILE: `pattern.R`
test_pattern <- function(model="tintori") {
	source("utils/pattern.R")
	source("utils/feature_selection.R")
	source("utils/utils.R")
	require(caret)
	require(MAST)
	require(scater)
	options(expressions = 5e5, warn=-1)
	verboseIt("-- START TEST")
	verboseIt("TEST binary logistic regression")
	data(infert)
	tl = 1:nrow(infert) %in% sample(1:nrow(infert), 1/3*nrow(infert))
	infert$induced <- as.numeric(infert$induced > 0)
	trainingSet <- infert[tl,]
	validationSet <- infert[!tl,]
	verboseIt("Single test")
	res <- logistic_regression(list(levels=trainingSet$case, 
			values=trainingSet$induced), 
		validationSet=list(levels=validationSet$case, 
			values=validationSet$induced), print=T)
	verboseIt("Cross validation")
	res <- cross_validation(k=8, infert$case, as.numeric(infert$induced > 0), print=T)
	verboseIt("TEST find gene expression patterns")
	if (model=="ciona") {
		load("../datasets/clust_ciona.Rdata")
		## remove housekeeping genes                 ##
		sceset <- ciona[which(!grepl("^ERCC-", rownames(ciona))), ]
	}
	if (model == "tintori") {
		load("../datasets/clust_tintori.Rdata")
		sceset <- tintori
	}
	if (model=="newciona") {
		load("../datasets/clust_newciona.Rdata")
		## remove housekeeping genes                 ##
		sceset <- cionan[which(!grepl("^ERCC-", rownames(cionan))), ]
	}
	## Data pre-processing        ##
	sceset <- updateSCESet(sceset)
	sceset <- sceset[sample(1:nrow(sceset), 2e3), sample(1:ncol(sceset), 15)]
	replicates <- NULL
	M <- filterMatrix(counts(sceset))
	M <- normalizeLogData(M)
	sca <- matrix2sca(M)
	sca <- filterGE(sca, X=6)
	sca <- featureSelection(sca, thres=0.75)
	#corMatrix <- cor(t(assay(sca)), method="spearman")
	#sca <- sca[!(1:nrow(sca) %in% findCorrelation(corMatrix, cutoff=0.75)), ]
	sca_r <- group_by_condition(sca, replicates)
	res <- getPatterns(sca_r, nbinsden=1e2, fast=F)
	mat <- as.matrix(dist(t(res$pattern), method="euclidean")/nrow(res$pattern))
	colnames(mat) <- sapply(colnames(mat), function(e) strsplit(e, "_")[[1]][2])
	rownames(mat) <- colnames(mat)
	heatmap_c(mat)
	verboseIt("-- END TEST")
}

#____________________#
# Clustering         #
#____________________#

#' Print clustering result with PCA
#'
#' Prints clustering result with PCA: clustering as colour
#' true labelling as point labels
#'
#' @param labels clustering
#' @param res a R named list with `matrix`: data matrix on which 
#'		the clustering has been performed
#' @param truth optional true labelling
#' @return g ggbiplot object
#'
#' @importFrom ggbiplot ggbiplot
#' @importFrom stats prcomp
#' @importFrom stats var
#'
#' @export
printCLUST <- function(labels, res, truth=NULL) {
	require(ggbiplot)
	## Performs PCA on data matrix ##
	pcaa <- prcomp(t(res$matrix[!apply(res$matrix, 1, function(e) var(e)==0), ]), center=T, scale=T)
	## Builds plot                 ##
	g <- ggbiplot(pcaa, obs.scale = 1, var.scale = 1, var.axes = F, 
	      groups = sapply(labels, function(e) paste0("Cluster #", e)), ellipse = TRUE, 
	      labels = truth, labels.size = 3,
	      circle = TRUE)
	print(g)
}

#' Tests consensus of correlation clusterings (-modified- Ailon et al.'s algorithm)
#'
#' @param n number of clusters
#' @param d dimension
#' 
#' FUNCTION: `clusterize`
#' FILE: `utils.R`
test_clustering <- function(n, d, ailon=F, times=1, percentile=0.15) {
	source("utils/utils.R")
	require(mclust)
	f <- function(n, d) return(do.call(cbind, unlist(lapply(1:n, function(i) 
		lapply(1:sample(2:3, 1), function(m) rnorm(d, mean=i, sd=i*1e-1))), recursive=F)))
	M <- f(n, d)
	truth <- apply(M, 2, function(x) round(mean(x)))/10
	applyClust <- function() return(clusterize(as.matrix(dist(t(M), method="euclidean"))[,], ailon=ailon, percentile=percentile)$clusters)
	print("KMEANS")
	print(system.time(kmeansCL <- lapply(1:times, function(i) kmeans(t(M), centers=n, iter.max=times)$cluster))[3])
	print("AILON")
	print(system.time(ailonCL <- lapply(1:times, function(i) applyClust()))[3])
	arisKMEANS <- sapply(1:times, function(i) adjustedRandIndex(kmeansCL[[i]], truth))
	cat(paste0("Median ARI value with ground truth over all ", times, " iterations: ", 
		round(median(arisKMEANS), 2), "\n"), "", file="")
	cat(paste0("Standard deviation of ARI value: ", 
		round(sqrt(var(arisKMEANS), 2)), "\n"), "", file="")
	print("AILON")
	arisAILON <- sapply(1:times, function(i) adjustedRandIndex(ailonCL[[i]], truth))
	cat(paste0("Median ARI value with ground truth over all ", times, " iterations: ", 
		round(median(arisAILON), 2), "\n"), "", file="")
	cat(paste0("Standard deviation of ARI value: ", 
		round(sqrt(var(arisAILON), 2)), "\n"), "", file="")
	cat(paste0("Mean cluster number difference (#AILON-#TRUTH): ", 
		mean(sapply(1:times, function(i) length(table(ailonCL[[i]]))-length(table(truth)))), "\n"), "", file="")
	printCLUST(ailonCL[[which.max(arisAILON)[1]]], list(matrix=M), truth)
	return(NULL)
}

###################################
# BONUSES                         #
###################################

#_____________________________________________#
# General tests                               #
#_____________________________________________#

test_triangular_inequality <- function(distMatrix) {
	n <- nrow(distMatrix)
	s <- apply(expand.grid(1:n, 1:n, 1:n), 1, function(x) {
		distMatrix[x[1], x[2]]+distMatrix[x[2], x[3]]>distMatrix[x[1], x[3]]
		})
	print(all(s))
	return(s)
}

#_____________________________________________#
# Tests for gene expression distribution      #
#_____________________________________________#

## Using psych::pairs.panels : 1 < n, m < 11

processMatrix <- function(sceset, dropout) {
	library(psych)
	mat <- counts(updateSCESet(sceset))
	return(switch(dropout,
			"1" = mat[rowSums(mat) == 0, ],
			"2" = mat[rowSums(mat) > 0, ],
			default = mat))
}

## Plots dependence b/w n randomly selected genes in 1 cell         ##
marginalDistr <- function(sceset, n=5, dropout="2") {
	M <- t(processMatrix(sceset, dropout))
	return(pairs.panels(M[sample(1:nrow(M), 2), sample(1:ncol(M), n)]))
}

## Plots dependence b/w m randomly selected cells for n randomly    ##
## selected                                                         ##
copulaDistr <- function(sceset, n=5, m=5, dropout="2") {
	M <- t(processMatrix(sceset, dropout))
	return(pairs.panels(M[sample(1:nrow(M), n), sample(1:ncol(M), m)]))
}

getList <- function(field, param, interval, sceset) {
	tmp <- unlist(lapply(interval, function(i) if (pData(sceset)[i, field] == param) i))
	return(tmp[!(is.null(tmp))])
}

## Plots dependence b/w m randomly selected cells for n randomly    ##
## selected respect to one parameter: embryo, cell ID, stage...     ##
copulaDistrSelected <- function(sceset, n=5, m=5, em=NULL, ct=NULL, stage=NULL, dropout="2") {
	tmp <- 1:ncol(sceset)
	if (!is.null(em)) tmp <- getList("Embryo", em, tmp, sceset)
	if (!is.null(ct)) tmp <- getList("Cell.ID", ct, tmp, sceset)
	if (!is.null(stage)) tmp <- getList("EmbryoStage", stage, tmp, sceset)
	return(copulaDistr(sceset[, tmp], n, m, dropout))
}

## Plots m-cell dependence for 1 randomly selected gene             ##
marginalDistrAll <- function(sceset, type="log2(TPM+1)", m=5, dropout="2") {
	M <- t(processMatrix(sceset, dropout))
	gene <- sample(1:nrow(M), 1)
	plot(1:ncol(M), M[gene, ], 
		main=sprintf("Gene expression in the dataset for gene %s", rownames(M)[gene]),
		xlab = "Cells", ylab = sprintf("%s counts", type), xaxt='n')
		axis(side=1, at=1:ncol(M), labels=colnames(M))
	return(NULL)
}

#_____________________________________________#
# Tests for binary logistic regression        #
#_____________________________________________#

sceset2sca <- function(sceset) {
	library(splatter)
	source("mineClusters.R")
	sca <- matrix2sca(counts(updateSCESet(sceset)))
	return(sca)
}

getFitData <- function(sceset, event="dropout") {
	library(MAST)
	sca <- sceset2sca(sceset)
	p <- nrow(sca)
	m <- ncol(sca)
	capture.output(res <- thresholdSCRNACountMatrix(assay(sca), 
				nbins=floor(nrow(sca)/1e2), min_per_bin=2))
	assay(sca) <- res$counts_threshold
	switch(event,
		"dropout" = {
			values <- rowSums(apply(assay(sca), 2, function(x) 
					if (sum(x) == 0) rep(0, p) else x/sum(x)))
			levels <- rep("", p)
			idx <- values == 0
			levels[idx] <- "dropout"
			levels[!idx] <- "not_dropout"
		},
		"mild" = {
			threshold <- sapply(as.character(res$bin), function(bin) res$cutpoint[bin])
			values <- assay(sca)
			levels <- apply(values, 2, function(x) ifelse(x > threshold, "high", "mild"))
		},
		default = {values <- NULL; levels <- NULL}
	)
	if (!is.null(levels) && !is.null(values)) {
		levels <- ifelse(levels == levels[1], 1, 0)
		values <- values-mean(values)
	}
	return(list(levels=levels, values=values))
}

## Tests linearity between dependent variables (dropout or expressiveness) and logits      ##
##   levels = 1/(1+exp(a*centered_values)) (levels being 0 or 1 vector)                    ##
testLinearity <- function(sceset, event=NULL, res=NULL) {
	sca <- sceset2sca(sceset)
	if (is.null(res) && !is.null(event)) res <- getFitData(sceset, event)
	levels <- res$levels
	values <- res$values
	linearModel <- glm(levels ~ values, family="binomial")
	nullModel <- glm(levels ~ 1, family="binomial")
	## McFadden R^2 index ~ R^2 coefficient for linear regression     ##
	## 1-logLik(linearModel)/logLik(nullModel)                        ##
	## = pscl::pR2(linearModel)[["McFadden"]]                         ##
	## Should be close to 1 for a good fit                            ##
	res <- list(Coefs=linearModel$coefficients, 
			Pvalue=1-logLik(linearModel)/logLik(nullModel), 
			Anova=anova(linearModel, test="Chisq"))
	plot(values, levels, col="green", 
		type="p",
		main=sprintf("Plotting centered gene expression values against %s events", event),
		xlab="centered gene expression values", 
		ylab=sprintf("%s events", event),
		ylim=c(0, 1),
		xlim=c(min(values), quantile(values, prob=0.25)))
	abline(res[[1]][1], res[[1]][2], col="red")
	return(res)
}

## Tests independence between residuals                                                    ##
testIndependence <- function(sceset, event=NULL, res=NULL) {
	library(lmtest)
	sca <- sceset2sca(sceset)
	if (is.null(res) && !is.null(event)) res <- getFitData(sceset, event)
	y <- res$levels
	x <- res$values
	linearModel <- glm(y ~ x, family="binomial")
	plot(linearModel, which=1)
	residuals <- linearModel$residuals[-length(linearModel$residuals)]
	last_residual <- linearModel$residuals[-1]
	modelres = lm(residuals ~ last_residual)
	residuals <- linearModel$residuals
	return(list(Test=dwtest(y ~ x), 
		Model=linearModel, 
		CHISQ=chisq.test(abs(linearModel$residuals)), 
		ResidualModel=modelres, 
		ACF=acf(residuals)
		))
}
