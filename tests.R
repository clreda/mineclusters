####################################################################
#                                                                  #
#               TESTS FILE                                         #
#                                                                  #
####################################################################

########################
## UTILS              ##
########################

#' Print messages
#'
#' Prints messages during the tests
#'
#' @param msg a R character string
#' @return msg
#'
#' @export
verboseIt <- function(msg) print(paste0("MSG: ", msg, "!"))

#' Get data from a SCESet object
#'
#' Gets matrix from a SCESet object
#'
#' @param sceset a R object of type SCESet
#' @return matrix data matrix from the SCESet
#'
#' @importFrom scater counts
#' @importFrom scater fpkm
#' @importFrom scater tpm
#' @importFrom scater cpm
#' 
#' @export
getData <- function(sceset) {
	require(scater)
	cts <- tryCatch({counts(sceset)}, error=function(e) NULL)
	fts <- tryCatch({fpkm(sceset)}, error=function(e) NULL)
	tts <- tryCatch({tpm(sceset)}, error=function(e) NULL)
	cpts <- tryCatch({cpm(sceset)}, error=function(e) NULL)
	value <- c(cts, fts, tts, cpts)
	value <- matrix(value, ncol=ncol(sceset), nrow=nrow(sceset))
	colnames(value) <- colnames(sceset)
	rownames(value) <- rownames(sceset)
	return(value)
}

#' Load model from Rdata file
#' 
#' Loads SCESet model from Rdata file
#'
#' @param model a R character string for the filename
#' @param path a R character string for the absolute/relative path of the file
#' @return sceset a SCESet object
#' 
#' @importFrom scater pData
#'
#' @export
loadModel <- function(model, path="../datasets/") {
	require(scater)
	## Removes housekeeping genes             ##
	removeERCC <- function(sceset) return(sceset[!grepl("^ERCC-", rownames(sceset)), ])
	load(paste0(path, "clust_", model, ".Rdata"))
	if (model == "newciona") model <- "cionan"
	## The SCESet object stored in the file   ##
	## should be named `model`                ##
	sceset <- get(model)
	if (grepl("ciona", model)) {
		sceset <- removeERCC(sceset)
	}
	if (model == "cionan") {
		## Creates `Cell.ID` and `Embryo` slots   ##
		pData(sceset)$Cell.ID <- sapply(pData(sceset)$cell, function(x) strsplit(as.character(x), "_")[[1]][2])
		pData(sceset)$Embryo <- sapply(pData(sceset)$cell, function(x) strsplit(as.character(x), "_")[[1]][1])	
	}
	cat(paste0("Dimensions ", model, " SCESet = ", nrow(sceset), " genes x ", ncol(sceset), " samples"), "\n", file="")
	tryCatch({return(updateSCESet(sceset))}, error=function(e) return(sceset))
}

#' Get a subset of the SCESet and the associated ground truth 
#'
#' Gets a subset of the SCESet and the associated ground truth 
#'
#' @param model a string character: name of the model stored in directory located at path
#' @param nr number of randomly selected lines in the final SCESet
#' @param nc number of randomly selected columns in the final SCESet
#' @param cutoff if not NULL the correlation threshold for the feature selection
#' @param sceset an optional SCESet object if the desired SCESet has already been created
#' @param path a R character string: absolute or relative path to the Rdata file named `model.Rdata`
#' @return res a R named list `sceset`: created, trimmed, no-zero-variance-line SCESet from the SCESet of model
#' 		`M`: the associated data matrix
#' 		`replicates`: optional labelling of the samples according to their biological/technical replicate number
#'		`vcell`: the ground truth to clustering
#' 		`vembryo`: another labelling to clustering
#'
#' @importFrom stats var
#' @importFrom stats cor 
#' @importFrom caret findCorrelation
#'
#' @export
getValues <- function(model, nr=NULL, nc=NULL, cutoff=NULL, sceset=NULL, path="../datasets/") {
	require(caret)
	if (is.null(sceset)) {
		## Loads SCESet                ##
		sceset <- loadModel(model, path=path)
		## Selects features            ##
		features <- if (!is.null(nr)) {
			sampleSet <- (1:nrow(sceset))[apply(getData(sceset), 1, function(x) var(x) != 0)]
			sample(sampleSet, min(nr, length(sampleSet)))
			} else 1:nrow(sceset)
		## Selects conditions          ##
		conditions <- if (!is.null(nc)) {
			conditions <- sample(1:ncol(sceset), min(nc, ncol(sceset)))
			} else 1:ncol(sceset)
		## Trims SCESet accordingly    ##
		sceset <- sceset[features, conditions]
		## Performs feature selection  ##
		if (!is.null(cutoff)) {
			corMatrix <- cor(getData(sceset), method="spearman")
			features <- (1:nrow(sceset))[!(1:nrow(sceset) %in% findCorrelation(corMatrix, cutoff=cutoff))]
			sceset <- sceset[features, ]
		}
	}
	M <- getData(sceset)
	conditions <- unique(sapply(colnames(sceset), function(e) strsplit(e, "[.]")[[1]][1]))
	names(conditions) <- 1:length(conditions)
	## Gets replicate numbers ##
	replicates <- as.vector(sapply(colnames(sceset), 
		function(e) names(conditions)[grepl(strsplit(e, "[.]")[[1]][1], conditions)]))
	## Gets true labelling    ##
	if (!is.null(model) && (grepl("ciona", model)||model=="tintori")) {
		vcell <- sceset$Cell.ID
		vembryo <- sceset$Embryo
	}
	else {
		vcell <- sceset$ref_cl
		vembryo <- if (all(sceset$cell_type2==vcell)) sceset$cell_type1 else sceset$cell_type2
	}
	return(list(sceset=sceset, M=M, replicates=replicates, vcell=vcell, vembryo=vembryo))
}

########################
## PLOTTING           ##
########################

#' Plot time graph for multiple runs of MineClusters
#'
#' Plots time graph for multiple runs of MineClusters with different dataset sizes
#' (in red) and a linear time line according to the runtime for the first test (in green)
#'
#' @param res the output of function `test_time`
#' 
#' @export
plot_time <- function(res) {
	## Slope of the line in green  ##
	slope=res$timers[1]/res$nr[1]
	## Plots time graph            ##
	plot(if (length(res$nc)==1) res$nr else res$nc, res$timers, type="l", col="red", 
		main=paste0("Runtime for different sizes in model ", res$model),
		xlim=if (length(res$nc==1)) c(min(res$nr), max(res$nr)) else c(min(res$nc), max(res$nc)),
		xlab=ifelse(length(res$nc)==1, "e3 genes", "e1 samples"), ylab="Time in sec.")
	## Plots linear time line      ##
	lines(if (length(res$nc)==1) res$nr else res$nc, slope*res$nr, col="green")
}

#' Plot heatmap for dissimilarity matrix
#'
#' Plots heatmap for dissimilarity matrix 
#' (can be ordered according to the ground truth)
#'
#' @param dissMatrix the dissimilarity matrix
#' @param truth the vector of dissimilarity matrix rows/columns 
#' 	corresponding to ground truth
#' @return heatmap heatmap object
#' 
#' @importFrom NMF aheatmap
#' 
#' @export
aheatmap_c <- function(dissMatrix, truth=NULL) {
	require(NMF)
	## Optional ordering of matrix     ##
	if (!is.null(truth)) rows <- rownames(dissMatrix)[order(truth)]
	else rows <- 1:nrow(dissMatrix)
	ordered <- ifelse(!is.null(truth), NA, T)
	dissMatrix <- dissMatrix[rows, rows]
	if (!is.null(truth)) {
		colnames(dissMatrix) <- sort(truth)
		rownames(dissMatrix) <- sort(truth)
	}
	return(aheatmap(dissMatrix, breaks=median(dissMatrix), border_color = NA, Rowv = ordered, Colv = ordered))
}

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
	pcaa <- prcomp(t(res$matrix[!apply(res$matrix, 1, function(e) var(e)==0), ]), center=TRUE, scale=TRUE)
	## Builds plot                 ##
	g <- ggbiplot(pcaa, obs.scale = 1, var.scale = 1, var.axes = F, 
	      groups = sapply(labels, function(e) paste0("Cluster #", e)), ellipse = TRUE, 
	      labels = truth, labels.size = 3,
	      circle = TRUE)
	print(g)
}

########################
## BENCHMARKS         ##
########################

#' Quick benchmark of mineClusters
#'
#' Quick benchmark of mineClusters on a subset of a given SCESet with CMEANS
#' KMEANS HCLUST DIANA SC3
#'
#' @param model a R character string
#' @param nr if not null the number of randomly selected features
#' @param nc if not null the number of randomly selected conditions
#' @param thres the correlation threshold for clustering (or a cluster number)
#' @param sceset an optional SCESet object on which the benchmark should be performed
#' @param thres_fs the feature correlation threshold for feature selection in mineClusters
#' @param fast logical for fast (but less accurate) version of mineClusters using only patterns
#' @param parallel logical for parallelization in mineClusters
#' @param path relative path to the imported datasets
#' @param use_gene_filter logical for gene filtering in mineClusters
#' @return res a R named list `dissimilarity`: dissimilarity matrix computed by mineClusters
#'				`pattern`: pattern matrix computed in mineClusters
#' 				`sceset`: SCESet object on which clustering has been performed
#'				`dendrogram`: HCLUST object returned by mineClusters
#' 				`comp`: named list of `clusterings`: named list of clusterings for each algorithm
#'							`statistics`: data frame of ARI results and runtimes for each algorithm
#'
#' @importFrom scater rowData
#' @importFrom scater logcounts
#' @importFrom mclust adjustedRandIndex
#' @importFrom e1071 cmeans
#' @importFrom fastcluster hclust
#' @importFrom cluster diana
#' @importFrom SC3 sc3
#'
#' @export
test_main <- function(model, nr=NULL, nc=NULL, thres=NULL, sceset=NULL, thres_fs=0.75, X=6, fast=F, parallel=T, path="../datasets/", use_gene_filter=T) {
	source("mineClusters.R")
	require(scater)
	require(mclust)
	require(e1071)
	require(fastcluster)
	require(cluster)
	require(SC3)
	verboseIt("-- START TEST")
	## Builds test SCESet from Rdata file in path          ##
	res <- getValues(model=model, sceset=sceset, nr=nr, nc=nc, cutoff=NULL, path=path)
	## Is the model already normalized? (in RPKM, etc.)    ##
	normalizedmodel <- !(!is.null(model)&&(model=="tintori"||grepl("ciona", model)))
	verboseIt(paste(model, "model is", ifelse(normalizedmodel, "already", "not"), "normalized"))
	## Gets matrix and SCESet                              ##
	M <- res$M
	sceset <- res$sceset
	rowData(sceset)$feature_symbol <- rownames(M)
	replicates <- NULL
	valuesCell <- res$vcell
	valuesEmbryo <- res$vembryo
	## Initializes `logcounts` slot                        ##
	logcounts(sceset) <- normalizeData(M)
	## True labelling                                      ##
	truth = table(valuesCell)[table(valuesCell)!=0]
	## Expected number of clusters                         ##
	nclusters = length(table(valuesCell)[table(valuesCell)!=0])
	## List of tested algorithms                           ##
	algos <- c("Rand", "HCLUST_avg", "HCLUST_cplt", "HCLUST_avg_cor", "HCLUST_cplt_cor", 
			"KMEANS", "CMEANS", "MINECLUSTERS", "DIANA", "SC3", "DIANA_cor")
	## Array of associated runtimes                        ##
	times <- rep(0, length(algos))
	## TEST: MINECLUSTERS                                  ##
	times[8] <- system.time(res <- mineClusters(M, options=list(thres=thres, replicates=replicates,
			 thres_fs=thres_fs, normalized=normalizedmodel, X=X, islog=F, fast=fast, parallel=parallel, use_gene_filter=use_gene_filter)))[3]
	if (is.null(res)) {
		verboseIt("Abort benchmark: critical error in mineClusters")
		return(NULL)
	}
	dissMatrix <- res$dissimilarity
	print(paste0("Expected cluster # = ", nclusters))
	print(paste0("Cluster # found = ", length(unique(res$clusters))))
	## Testing functions                                   ##
	fcts <- list(
		## RANDOM clustering                           ##
		function(ls, d) { ls[[1]] <- 1:ncol(M); return(ls) },
		## Average-linkage hierarchical clustering          ##
		## on dissimilarity matrix computed by mineClusters ##
		function(ls, d) { ls[[2]] <- cutree(fastcluster::hclust(d, method="average"), k=nclusters); return(ls) },
		## Complete-linkage hierarchical clustering         ##
		## on dissimilarity matrix computed by mineClusters ##
		function(ls, d) { ls[[3]] <- cutree(fastcluster::hclust(d, method="complete"), k=nclusters); return(ls) },
		## Average-linkage hierarchical clustering     ##
		## on distance matrix computed by correlation  ##
		function(ls, d) { dcorMatrix <- as.dist(1-cor(logcounts(sceset), method="spearman"))
			ls[[4]] <- cutree(fastcluster::hclust(d, method="average"), k=nclusters)
			return(ls) },
		## Complete-linkage hierarchical clustering    ##
		## on distance matrix computed by correlation  ##
		function(ls, d) { dcorMatrix <- as.dist(1-cor(logcounts(sceset), method="spearman"))
			ls[[5]] <- cutree(fastcluster::hclust(d, method="complete"), k=nclusters)
			return(ls) },
		## K-MEANS                                     ##
		function(ls, d) { ls[[6]] <- kmeans(t(logcounts(sceset)), iter.max=100, centers=nclusters)$cluster; return(ls) },
		## Fuzzy K-MEANS                               ##
		function(ls, d) { ls[[7]] <- cmeans(t(logcounts(sceset)), iter.max=100, centers=nclusters)$cluster; return(ls) },
		## MineClusters                                ##
		function(ls, d) { ls[[8]] <- res$clusters; return(ls) },
		## DIANA on dissimilarity matrix in mineClusters    ##
		function(ls, d) { ls[[9]] <- cutree(as.hclust(diana(d)), k=nclusters); return(ls) },
		## SC3                                         ##
		function(ls, d) { 
			ls[[10]] <- sc3(sceset, gene_filter = !normalizedmodel, 
					ks=nclusters)[[paste0("sc3_", nclusters, "_clusters")]]
			return(ls) },
		## DIANA on correlation matrix                 ##
		function(ls, d) { dcorMatrix <- as.dist(1-cor(logcounts(sceset), method="spearman"))
			ls[[11]] <- cutree(as.hclust(diana(dcorMatrix)), k=nclusters)
			return(ls) }
	)
	## Stores clusterings for each algorithm               ##
	ls <- rep(list(NULL), length(algos))
	names(fcts) <- algos
	names(ls) <- algos
	## Apply algorithms                                    ##
	for (i in (1:length(algos))) {
		cat(paste0("Algorithm #", i, ": ", algos[i]), "\n", file="")
		times[i] <- system.time(ls <- fcts[[i]](ls, as.dist(dissMatrix)))[3] + ifelse(i==8, times[8], 0)
		if (i %in% c(2, 3, 9, 11)) times[i] <- times[i] + times[8]
	}
	tests <- c("Cell.ID", "Other", "Time")
	## Stores statistics                                   ##
	ress <- matrix(0, ncol=length(algos), nrow=length(tests))
	rownames(ress) <- tests
	colnames(ress) <- algos
	if (length(ls) != length(algos)) {print(paste0("ERROR: ", length(algos), " != ", length(ls))); return(NULL)}
	for (i in 1:length(ls)) ress[, i] <- c(sapply(
			sapply(list(valuesCell, valuesEmbryo), function(d) 
					adjustedRandIndex(ls[[i]], d)), function(e) round(e, 3)), 
						times[i])
	print(as.data.frame(ress))
	## Ranking                                             ##
	for (i in 1:length(tests)) {
		cat(paste0("Ranking for ", rownames(ress)[i]), " = ")
		cat(paste0(colnames(ress)[sort(ress[i, ], 
			decreasing=(rownames(ress)[i]!="Time"), index.return=T)$ix], collapse=" > "), "\n")
	}
	res$matrix <- M
	verboseIt("-- END TEST")
	## Builds final result file                            ##
	res <- list(dissimilarity=res$dissimilarity, pattern=res$pattern, sceset=sceset, 
		dendrogram=res$dendrogram, comp=list(clusterings=ls, statistics=as.data.frame(ress)))
	save(res, file=paste0(c(Sys.time(), model, 
			paste0("nr=", ifelse(is.null(nr), "all", nr)), 
			paste0("nc=", ifelse(is.null(nc), "all", nc)),
			if (fast) "fast" else NULL,
			"benchmark.Rdata"), collapse="_"))
	vars <- ls(envir = globalenv())
 	rm(list = vars[vars != "res"], envir = globalenv())
	gc()
	cat("\014")
	source("tests.R")
	return(res)
}

#' Boxplot for multiple benchmarks
#' 
#' Boxplot for multiple benchmarks generated by `test_main`
#'
#' @param path a R character string for location of benchmark result files
#' @param results a R character string: name of an already built boxplot result file
#' @return res a boxplot result file that can be called again via `test_mean_perf`
#' 
#' @importFrom graphics boxplot
#'
#' @export
test_mean_perf <- function(path=".", results=NULL) {
	if (is.null(results)) {
		## Gets all benchmark result files created by `test_main`   ##
		datasets <- dir(path = path, pattern = '_benchmark.Rdata')
		## Aggregates statistics results                            ##
		results <- do.call(rbind, lapply(1:length(datasets), function(i) {
								load(file.path(path, datasets[i]))
								matrix(res$comp$statistics["Cell.ID", ], nrow=1)
							}))
		load(file.path(path, datasets[1]))
		colnames(results) <- colnames(res$comp$statistics)
		results2 <- do.call(rbind, lapply(1:ncol(results), function(i) {
								m <- matrix(results[, i], ncol=1)
							}))
		rn <- sapply(colnames(results), function(coli) rep(coli, nrow(results)))
		results2 <- cbind(results2, matrix(rn, ncol=1))
		colnames(results2) <- c("ARI", "Algorithms")
		results2 <- as.data.frame(results2)
		save(results2, file=paste0(Sys.time(), "_results.Rdata"))
	}
	else load(results)
	## Removes the results associated with the following                ##
	## algorithms in the boxplot                                        ##
	idx <- !(results2$Algorithms %in% c("Rand", "HCLUST_avg", "HCLUST_avg_cor", "HCLUST_cplt_cor", 
			"CMEANS", "DIANA_cor"))
	y <- unlist(results2$ARI)[idx]
	x <- unlist(results2$Algorithms)[idx]
	boxplot(y~x, main="Performance Boxplot (in ARI)", xlab="Algorithms", ylab="ARI value")
	gc()
	return(results2)	
}

#' Get runtimes for multiple runs of mineClusters on different subset sizes
#' of a given model
#'
#' Gets runtimes for multiple runs of mineClusters on different subset sizes
#' of a given model that can be then plotted with function `plot_time`
#'
#' @param model a R string character for the model to test
#' @param nr a vector of the feature number for each run (x 1e3)
#' @param nc a vector of the condition number for each run (x 10)
#' @param thres_fs the feature correlation threshold in mineClusters
#' @param parallel logical for parallelization
#' @param ret logical to save each run output in a file
#' @param path a R character string: location of the model Rdata file
#' @return res object that can be plotted using `plot_time`
#'
#' @export
test_time <- function(model="tintori", nr=NULL, nc=NULL, thres_fs=0.75, parallel=T, ret=F, path="../datasets/") {
	source("mineClusters.R")
	verboseIt("-- START TEST")
	times <- c()
	## Function called for each run         ##
	startIt <- function(r, c, times) {
		r = if (is.null(r)) NULL else r*1e3
		c = if (is.null(c)) NULL else c*10
		print_r = if (is.null(r)) "all" else r
		print_c = if (is.null(c)) "all" else c
		verboseIt(paste0("#rows = ", print_r, " & #cols = ", print_c))
		res <- getValues(model=model, nr=r, nc=c, cutoff=NULL, path=path)
		times <- c(times, system.time(res <- mineClusters(res$M, 0.75, thres_fs=thres_fs, parallel=parallel))[3])
		if (ret) save(res, file=paste0(c(Sys.time(), model, "nr=", r, "nc=", c, 
			ifelse(parallel, "parallel", ""), "result.Rdata"), collapse="_"))
		cat(paste0("Total Runtime = ", round(times[length(times)], 3), " s. \n"), "", file="")
		return(times)
	}
	for (r in nr) {
		for (c in nc) times <- startIt(r, c, times)
		if (is.null(nc)) times <- startIt(r, NULL, times)
	}
	if (is.null(nr)) {
		for (c in nc) times <- startIt(NULL, c, times)
		if (is.null(nc)) times <- startIt(NULL, NULL, times)
	}
	verboseIt("-- END TEST")
	res <- list(nr=nr, nc=nc, timers=times, model=model)
	save(res, file=paste0(Sys.time(), "_", res$model, "_nr=", paste0(res$nr, collapse="-"), 
			"_nc=", paste0(res$nc, collapse="-"), ".Rdata"))
	gc()
	return(res)
}
