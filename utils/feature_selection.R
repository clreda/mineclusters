##################################################################
#                                                                #
#   FEATURE SELECTION FOR MineClusters ALGORITHM                 #
#                                                                #
##################################################################

######################################################
# Normalization                                      #
######################################################

#' log2 quasi-TPM+1 normalization of matrix
#'
#' Converts a matrix into a log2 quasi-TPM+1 normalized matrix
#' It does NOT normalize for gene length (we do not possess this information),
#' but it accounts for sequencing depth.
#'
#' @param M raw counts matrix
#' @return M normalized matrix
#'
#' @export
normalizeData <- function(M) {
	return(as.matrix(apply(M, 2, function(x) log(x/sum(x)*10**6+1)/log(2))))
}

#' log2+1 transformation of matrix
#'
#' Converts a matrix into a log2 +1 transformed matrix.
#'
#' @param M normalized matrix
#' @return M regularized matrix
#' 
#' @export
normalizeLogData <- function(M) {
	return(as.matrix(log(M+1)/log(2)))
}

######################################################
# Filtering                                          #
######################################################

#' Filter matrix values
#' 
#' Filters NA values and null variance lines.
#'
#' @param M matrix
#' @return M filtered matrix
#' 
#' @export
filterMatrix <- function(M) {
	M <- na.omit(M)[,]
	## Delete features of null variance ##
	M <- M[!apply(M, 1, function(x) var(x)==0), ]
	return(M)
}

#' Filter genes on gene expression levels
#' 
#' Filters genes based on gene expression levels.
#'
#' @param sca SingleCellAssay object
#' @param X percentage of samples REF: SC3 article (Kiselev et al., 2017)
#'          (removes ubiquitous features: present in >= (100-X)% samples)
#'          (removes rare features: present in <= X% samples            )
#'          "Present" means here having a non-zero expression value, which should not
#'          rely on gene length (thus normalization)
#' @param lower 
#' @return sca filtered SingleCellAssay object
#'
#' @importFrom MAST freq
#' @importFrom SummarizedExperiment assay
#' 
#' @export
filterGE <- function(sca, X=6) {
	sca <- sca[which(freq(sca) > X*1e-2), ]
	sca <- sca[which(freq(sca) < 1-X*1e-2), ]
	return(sca)
}

######################################################
# Feature Selection based on Correlation             #
######################################################

#' Feature selection based on correlation
#'
#' Applies custom feature selection to remove genes that
#' are too positively correlated to others.
#'
#' @param sca SingleCellAssay object
#' @param thres threshold for similarity between features
#' @param use_caret logical for using filter method from caret package
#' @return sca trimmed SingleCellAssay object
#' 
#' @importFrom stats cor
#' @importFrom stats mean
#' @importFrom stats var
#' @importFrom stats dist
#' @importFrom SummarizedExperiment assay
#' @importFrom fastcluster hclust
#' @importFrom stats cutree
#' @importFrom stats findCorrelation
#' 
#' @export
featureSelection <- function(sca, thres, use_caret=(nrow(sca)>1e4)) {
	## Correlation with caret package                                       ##
	if (use_caret) {
		highlyCorrelated <- findCorrelation(cor(t(assay(sca)), method="spearman"), cutoff=thres, exact=F)
		return(sca[-highlyCorrelated, ])
	}
	## Custom method                                                        ##
	n <- nrow(assay(sca))
	means <- apply(assay(sca), 1, mean)
	stdevs <- apply(assay(sca), 1, function(x) sqrt(var(x)))
	upper <- means+stdevs
	lower <- means-stdevs
	## Genes are considered truly redundant iff. they are highly positively ##
	## correlated and their mean expression value is roughly similar        ##
	## redundant is a matrix of size n such as                              ##
	## redundant[i, j] = 1 iff. genes i and j have similar mean expression  ##
	## value i.e. |mean(GE_i)-mean(GE_j)| <= sd(GE_j)                       ##
	redundant <- apply(as.matrix(means, ncol=1), 1, function(x) (x >= lower)*(x <= upper))
	## Maximum distance value: 1                                            ##
	cordist <- (1-cor(t(assay(sca)), method="spearman"))/2
	## Non redundant genes have maximum distance value                      ##
	cordist[redundant==0] <- 2
	hc <- cutree(hclust(as.dist(cordist), method="complete"), h=1-thres)
	## Selection of the "first" elements of each cluster                    ##
	return(sca[match(1:max(hc), hc), ])
}

#################################
## ONE LINER FUNCTION          ##
#################################

preprocessing_step <- function(M, options) {
	clock(paste0("Initial dataset ", nrow(M), "x", ncol(M)))
	M <- filterMatrix(M)
	clock(paste0("Filtering... ", nrow(M), "x", ncol(M)))
	if (!options$islog) M <- if (options$normalized) normalizeLogData(M) else normalizeData(M)
	else M <- if (options$normalized) M else normalizeData(exp(M))
	clock(paste0("Normalization... ", nrow(M), "x", ncol(M)))
	sca <- matrix2sca(M)
	sca <- filterGE(sca, X=options$X)
	clock(paste0("Filtering GE... ", nrow(sca), "x", ncol(sca)))
	sca <- featureSelection(sca, thres=options$thres_fs)
	clock(paste0("Feature Selection... ", nrow(sca), "x", ncol(sca)))
	sca_replicates <- group_by_condition(sca, options$replicates)
	return(list(sca=sca, sca_replicates=sca_replicates))
}
