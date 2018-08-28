################################################################
#                                                              #
#               MineClusters ALGORITHM CORE CODE               #
#                                                              #
################################################################

source("import.R")

#' Cell clustering of gene expression matrix
#'
#' Clusters cell according to gene expression matrix.
#'
#' @param M (normalized) gene expression matrix
#' @param options R named list of
#' `thres` correlation threshold for clustering (if < 1) or number of clusters (otherwise)
#' `replicates` vector of size ncol(M) denoting experiment replicates for each condition
#'                   by anything that is convertible to an integer
#' `thres_fs` correlation threshold for feature selection
#' `X` frequency-based gene trimming criterion in sample percentage
#' `normalized` is matrix M already normalized? default FALSE
#' `islog` is matrix M already in log counts? default FALSE
#' `fast` logical for using a faster pipeline which does not compute the confidence coefficient
#'             using copulas
#' `parallel` logical for using multiple cores
#' `retParams` logical for returning copula and copula parameters (can be memory-consuming)
#' `use_gene_filter` logical for filtering genes
#' @return res a named R list that contains 
#'             "clusters": vector of labels for each sample
#'             "matrix": the trimmed normalized gene expression matrix
#'             "copulas": the list of copulas for each condition (group of samples) = list of cumulative prob. functions
#'             "dissimilarity": dissimilarity matrix computed with the copulas
#' 
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel clusterExport
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#'
#' @export
mineClusters <- function(M, options=NULL) {
    #__________________________________________________________________________________________________________________________#
    #                           PRE-PROCESSING                                                                                 #
    #__________________________________________________________________________________________________________________________#
	default <- list(thres=NULL, replicates=NULL, thres_fs=0.75, X=6, normalized=F, islog=F, fast=F, parallel=T, retParams=F, use_gene_filter=T)
	if (is.null(options)) options <- default
	else for (i in names(default)) if (is.null(options[[i]])) options[[i]] <- default[[i]]
	if (!is.null(options$replicates)) options$replicates <- normalize_labels(options$replicates)
	else options$replicates <- 1:ncol(M)
	if (options$use_gene_filter) timer <- system.time(res <- preprocessing_step(M, options))[3]
	else timer <- system.time({
		sca <- matrix2sca(M)
		res <- list(sca=sca, sca_replicates=group_by_condition(sca, options$replicates))
	})[3]
	if (is.null(res)) return(NULL)
	sca <- res$sca
	sca_replicates <- res$sca_replicates
	rm(list="res")
	gc()
	if (!options$parallel) ncores <- 1
	clock(paste0("Start with ", ncores, "/", nacores, " cores", ifelse(ncores==1, ": no multi-core", "")))
	clock(paste0("Dimension reduction ", nrow(sca), " x ", ncol(sca)))
	if (nrow(sca) == 0) {
		cat("ERROR: No feature has been selected. You should maybe increase input argument thres_fs.")
		return(NULL)
	}
    	clock(paste0("(1/", ifelse(options$fast, "3", "6"), ") Normalization/Feature Selection/Replicate Summary"), timer)

    #__________________________________________________________________________________________________________________________#
    #                           FIND MOST PROBABLE PATTERNS                                                                    #
    # <=> Classification for each cell of every gene into dropout, mildly-expressed and highly-expressed                       #
    #__________________________________________________________________________________________________________________________#

	timer <- system.time(res <- getPatterns(sca_replicates, options))[3]
	if (is.null(res)) return(NULL)
	model_dropout <- res$dropout
	model_mild <- res$mild
	P <- res$pattern
	Q <- res$threshold
	peaks <- res$peaks
	rm(list="res")
	gc()
    	clock(paste0("(2/", ifelse(options$fast, "3", "6"), ") Dropout Model/Expressiveness Model/Pattern Matrix"), timer)

	if (options$fast) {
		timer <- system.time(res <- fast_clustering(M, P, options))[3]
		clock("(3/3) Fast Clustering (w/o confidence coefficient)", timer)
		return(res)
	}
    
    #__________________________________________________________________________________________________________________________#
    #                       COMPUTATION OF THE PROBABILITY DISTRIBUTION                                                        #
    # <=> Computation of parameters for marginal and for bimodal Gaussian copula for each condition                            #
    #__________________________________________________________________________________________________________________________# 

	timer <- system.time(params <- parameter_estimation(sca, model_dropout, model_mild, peaks, P, options))[3]
	clock("(3/6) Parameter Estimation", timer)
        timer <- system.time({
		P <- P[, options$replicates]
		Q <- Q[, options$replicates]
		QQ <- apply(Q, 2, function(x) transformQ(x, params))
	})[3]
	clock("(4/6) Bimodal Copula Estimation", timer)
    
    #__________________________________________________________________________________________________________________________#
    #                       PATTERN MERGING                                                                                    #
    # <=> Compute distance value between each pair of cells                                                                    #
    #__________________________________________________________________________________________________________________________#

        timer <- system.time({
		m <- ncol(sca)
		dissMatrix <- FBM(m, m)
		diag(dissMatrix) <- 0
		if (ncores > 1) {
			tryCatch({
			cl <- makeCluster(ncores, type=c("PSOCK", "FORK")[1], outfile="")
			registerDoParallel(cl, ncores)
			clusterEvalQ(cl=cl, library(mnormt, warn.conflicts=F, quietly=T))
			clusterEvalQ(cl=cl, options(warn=-1))
			clusterExport(cl=cl, 
				varlist=c("stopCluster", "pmnormCustom", "sadmvnCustom", "ncores", "copula", "params"),
				envir=.GlobalEnv)
			tmp <- foreach(k = 1:(m-1)) %dorng% {
				cat(ifelse(k%%5==0, k, "."), "", file="")
				dissMatrix[k, (k+1):m] <- as.vector(apply(as.matrix((k+1):m, ncol=1), 1, function(l) {
					norm(P[, k]-P[ ,l], "2")**2/nrow(P)*(1-copula(QQ[, l], params, k)*copula(QQ[, k], params, l))
				}))
				gc()
				NULL
			}
			stopCluster(cl)
			},
			error=function(e) {
				if (e[[1]]=="Error in sadmvnCustom(lower = rep(-Inf, d), upper = x[j, ], mean[j, ],  :") {
					clock("ERROR: Too few selected features. Perhaps you should increase the value of thres_fs.\n")
				}
				else print(e[[1]])
				gc()
				stopCluster(cl)
				return(NULL)
			})
		}
		else {
			for (k in 1:(m-1)) {
				cat(ifelse(k%%5==0, k, "."), "", file="")
				dissMatrix[k, (k+1):m] <- as.vector(apply(as.matrix((k+1):m, ncol=1), 1, function(l) {
					norm(P[, k]-P[ ,l], "2")**2/nrow(P)*(1-copula(QQ[, l], params, k)*copula(QQ[, k], params, l))
				}))
				gc()
				NULL
			}
		}
		dissMatrix <- as.matrix(dissMatrix[,], nrow=m, ncol=m, byrow=T)
		dissMatrix <- dissMatrix + t(dissMatrix)
		colnames(dissMatrix) <- colnames(sca)
		rownames(dissMatrix) <- colnames(sca)
	})[3]
	cat("\n", "", file="")
	clock("(5/6) Distance Matrix Computation", timer)
    
    #__________________________________________________________________________________________________________________________#
    #                       CLUSTERING                                                                                         #
    #__________________________________________________________________________________________________________________________# 

	timer <- system.time(res <- regular_clustering(dissMatrix, sca, options, params, P))[3]
	clock("(6/6) Clustering", timer)
	return(res)
}
