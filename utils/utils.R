##################################################################
#                                                                #
#       GENERAL UTILS FOR MineClusters ALGORITHM                 #
#                                                                #
##################################################################

#' Print current runtime
#'
#' Prints the runtime.
#'
#' @param msg message to print with the runtime
#' @param timer elapsed time
#' 
#' @export
clock <- function(msg, timer=NULL) {
  if (!is.null(timer)) cat(sprintf("\n[TIME %s: %2.2f seconds]", msg, timer), "\n")
  else cat(sprintf("\n[MSG %s]", msg), "\n")
}

#' Normalize labels
#'
#' Change a vector of labels into a vector of integers from 1 to k, where
#' k is the total number of distinct labels.
#'
#' @param labels label vector
#' @return res normalized label vector
#' 
#' @export
normalize_labels <- function(labels) {
	u <- unique(labels)
	names(u) <- seq_len(length(u))
	return(as.integer(sapply(labels, function(x) names(u[u==x]))))
}

#' Convert a matrix into a correct SingleCellAssay
#'
#' Converts a matrix into a correct SingleCellAssay.
#'
#' @param M normalized matrix
#' @return sca associated SingleCellAssay
#' 
#' @importFrom MAST FromMatrix
#' @importFrom SummarizedExperiment assay
#' 
#' @export
matrix2sca <- function(M) {
	sca <- FromMatrix(M, 
		## Column metadata          ##
		data.frame(wellKey=colnames(M), row.names=colnames(M)), 
		## Feature metadata         ##
		data.frame(primerid=row.names(M), row.names=row.names(M))
	)
	return(sca)
}

#' Summarize matrix according to replicates of a same condition
#'
#' Summarizes matrix according to biological/technical replicates of a same condition.
#'
#' @param sca SingleCellAssay object
#' @param replicates the vector that labels samples of a same condition
#' @param sfun the gene expression aggregating function
#' @return sca SingleCellAssay object with aggregated expression values
#' 
#' @importFrom stats mean
#' 
#' @export
group_by_condition <- function(sca, replicates, sfun=mean) {
	if (is.null(replicates)) return(sca)
	M <- do.call(cbind, lapply(unique(replicates), function(label) {
				if (length(which(replicates == label)) > 1) {
					apply(assay(sca)[, replicates == label], 1, sfun)
				}
				else assay(sca)[, replicates == label]
		}))
	row.names(M) <- row.names(assay(sca))
	colnames(M) <- unique(replicates)
	return(matrix2sca(M))
}

#' Correlation clustering with cluster number inference 
#'
#' Correlation clustering with cluster number inference
#' as described in Ailon et al.'s (2005).
#'
#' @param dissMatrix distance or dissimilarity matrix with coefficients between 0 and 1
#' @param thres element pairwise correlation threshold
#' @param ailon logical to use Ailon et al.'s method (and consensus clustering)
#' @param iter.max maximum number of iterations for Ailon et al.'s method
#' @return res list: `clusters`: result of clustering
#' 			`pivots`: list of selected pivots for reproducibility
#' 
#' @export
clusterize <- function(dissMatrix, thres=NULL, ailon=F, iter.max=1e2) {
	m <- ncol(dissMatrix)
	## Default threshold is one of the percentiles of dissimilarity matrix ##
	thres <- quantile(dissMatrix, prob=ifelse(is.null(thres), 0.75, thres))
	## Corresponding affinity matrix                                       ##
	## affinityMatrix <- ((-1)**(as.numeric(m>thres)))*abs(m-thres)        ##
	pivots <- NULL
	if (ailon) {
		clusters <- matrix(0, nrow=iter.max, ncol=m)
		for (i in 1:iter.max) {
			V <- seq_len(m)
			while (length(V)>1) {
				pivot <- sample(V, 1)
				pivots <- c(pivots, pivot)
				## If dissMatrix is a proper dissimilarity matrix ##
				## i.e. dissMatrix[x, x] = 0, it should be OK     ##
				C <- c(V[dissMatrix[pivot, V]<thres], pivot)
				clusters[i, C] <- pivot
				V <- which(clusters[i, ]==0)
			}
			clusters[i, V] = V
		}
		## Consensus clustering ##
		consensus <- sapply(1:m, function(j) sapply(1:m, function(i) sum(clusters[, i] == clusters[, j])))/iter.max
		clusters <- rep(0, m)
		for (j in 1:m) clusters[consensus[j, ]>=0.5] <- j
	}
	else {
		clusters <- rep(0, m)
		V <- seq_len(m)
		distPts <- apply(dissMatrix, 2, sum)
		while (length(V)>1) {
			pivot <- V[which.max(distPts[V])[1]]
			pivots <- c(pivots, pivot)
			## If dissMatrix is a proper dissimilarity matrix ##
			## i.e. dissMatrix[x, x] = 0, it should be OK     ##
			C <- c(V[dissMatrix[pivot, V]<thres], pivot)
			clusters[C] <- pivot
			V <- which(clusters==0)
		}
		clusters[V] = V
	}
	return(list(clusters=normalize_labels(clusters), pivots=pivots))
}

###########################################################################################
## Slightly faster versions of pre-existing functions in packages "mnormt" and "Matrix"  ##
###########################################################################################

upperTriCustom <- function(x) return(.Internal(row(dim(x))) < .Internal(col(dim(x))))

## mnormt::sadmvn
sadmvnCustom <- function (lower, upper, mean, varcov, utrho, maxpts = 2000*d, abseps = 1e-06, releps = 0) {
    upperTriCustom <- function(x) return(.Internal(row(dim(x))) < .Internal(col(dim(x))))
    if (any(lower > upper)) 
        stop("lower>upper integration limits")
    if (any(lower == upper)) 
        return(0)
    d <- as.integer(if (is.matrix(varcov)) ncol(varcov) else 1)
    sd <- rep(1, ncol(varcov))
    rho <- varcov
    lower <- as.double((lower - mean)/sd)
    upper <- as.double((upper - mean)/sd)
    if (d == 1) 
        return(pnorm(upper) - pnorm(lower))
    infin <- rep(2, d)
    infin <- replace(infin, (upper == Inf) & (lower > -Inf), 
        1)
    infin <- replace(infin, (upper < Inf) & (lower == -Inf), 
        0)
    infin <- replace(infin, (upper == Inf) & (lower == -Inf), 
        -1)
    infin <- as.integer(infin)
    if (any(infin == -1)) {
        if (all(infin == -1)) 
            return(1)
        k <- which(infin != -1)
        d <- length(k)
        lower <- lower[k]
        upper <- upper[k]
        if (d == 1) 
            return(pnorm(upper) - pnorm(lower))
        rho <- rho[k, k]
	utrho <- upperTriCustom(as.matrix(rho))
        infin <- infin[k]
        if (d == 2) 
            return(biv.nt.prob(0, lower, upper, rep(0, 2), rho))
    }
    lower <- replace(lower, lower == -Inf, 0)
    upper <- replace(upper, upper == Inf, 0)
    correl <- as.double(rho[utrho])
    maxpts <- as.integer(maxpts)
    abseps <- as.double(abseps)
    releps <- as.double(releps)
    error <- as.double(0)
    value <- as.double(0)
    inform <- as.integer(0)
    result <- .Fortran("sadmvn", d, lower, upper, infin, correl, 
        maxpts, abseps, releps, error, value, inform, PACKAGE = "mnormt")
    prob <- result[[10]]
    attr(prob, "error") <- result[[9]]
    attr(prob, "status") <- switch(1 + result[[11]], "normal completion", 
        "accuracy non achieved", "oversize")
    return(prob)
}

## mnormt::pmnorm
pmnormCustom <- function (x, mean = rep(0, d), varcov, utrho, ...) {
    d <- NCOL(varcov)
    x <- if (is.vector(x)) 
        matrix(x, 1, d)
    else data.matrix(x)
    n <- nrow(x)
    if (is.vector(mean)) 
        mean <- outer(rep(1, n), as.vector(matrix(mean, d)))
    if (d == 1) 
        p <- as.vector(pnorm(x, mean, sqrt(varcov)))
    else {
        pv <- numeric(n)
        for (j in 1:n) p <- pv[j] <- if (d == 2) 
            biv.nt.prob(0, lower = rep(-Inf, 2), upper = x[j, 
                ], mean[j, ], varcov)
        else sadmvnCustom(lower = rep(-Inf, d), upper = x[j, ], mean[j, 
            ], varcov, utrho, ...)
        if (n > 1) 
            p <- pv
    }
    return(p)
}

## Matrix::nearPD
nearPDCustom <- function (x, eig.tol = 1e-06, conv.tol = 1e-07, posd.tol = 1e-08, maxit = 100, 
    conv.norm.type = "I") {
    n <- ncol(x)
    D_S <- x
    D_S[] <- 0
    X <- x
    iter <- 0
    converged <- FALSE
    conv <- Inf
    while (iter < maxit && !converged) {
        Y <- X
        R <- Y - D_S
        e <- eigen(R, symmetric = TRUE)
        Q <- e$vectors
        d <- e$values
        p <- d > eig.tol * d[1]
        if (!any(p)) 
            stop("Matrix seems negative semi-definite")
        Q <- Q[, p, drop = FALSE]
        X <- tcrossprod(Q * rep(d[p], each = nrow(Q)), Q)
        D_S <- X - R
        diag(X) <- 1
        conv <- norm(Y - X, conv.norm.type)/norm(Y, conv.norm.type)
        iter <- iter + 1
        converged <- (conv <= conv.tol)
    }
    if (!converged) 
        warning(gettextf("'nearPD()' did not converge in %d iterations", 
            iter), domain = NA)
    e <- eigen(X, symmetric = TRUE)
    d <- e$values
    Eps <- posd.tol * abs(d[1])
    if (d[n] < Eps) {
            d[d < Eps] <- Eps
            Q <- e$vectors
            o.diag <- diag(X)
            X <- Q %*% (d * t(Q))
            D <- sqrt(pmax(Eps, o.diag)/diag(X))
            X[] <- D * X * rep(D, each = n)
    }
    diag(X) <- 1
    return(matrix(as.vector(X), n, n))
}

##########################################
## ONE LINER for fast clustering step   ##
##########################################

fast_clustering <- function(M, P, options) {
	dissMatrix <- as.matrix(dist(t(P), method="euclidean")/nrow(P))
	rownames(dissMatrix) <- colnames(M)
	colnames(dissMatrix) <- colnames(M)
	nclusters <- if (is.null(options$thres)||(options$thres < 1)) {
		length(table(clusterize(dissMatrix, thres=options$thres)$clusters))
	} else options$thres
	hc <- hclust(as.dist(dissMatrix), method="complete")
	clusters <- cutree(hc, k=nclusters)
	res <- list(clusters=clusters, 
			dendrogram=as.dendrogram(hc),
			matrix=M, 
			pattern=P,
			dissimilarity=dissMatrix)
	return(res)
}

##########################################
## ONE LINER for regular clustering step##
##########################################

regular_clustering <- function(dissMatrix, sca, options, params, P) {
	nclusters <- if (is.null(options$thres)||(options$thres < 1)) {
		length(table(clusterize(dissMatrix, thres=options$thres)$clusters)) 
	} else options$thres
	hc <- hclust(as.dist(dissMatrix), method="complete")
	clusters <- cutree(hc, k=nclusters)
	res <- c(list(clusters=clusters, dendrogram=as.dendrogram(hc)),
			list(matrix=assay(sca), 
			ifelse(options$retParams, 
				list(copulas=list(params=params)), 
				list(NULL)), 
			pattern=P,
			dissimilarity=dissMatrix))
	return(res)
}
