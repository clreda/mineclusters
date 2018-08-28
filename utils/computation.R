##################################################################
#                                                                #
#                COMPUTATIONS USED IN MineClusters               #
#                                                                #
##################################################################

sumSquares <- function(x) return(x %*% x)

######################################
# Dispersion parameters estimation   #
######################################

#' Estimate value of delta for Negative Binomial distr.
#'
#' Estimates the value of the dispersion parameter delta for Negative Binomial distr.
#'
#' @param M data matrix
#' @return delta a vector of float values
#' 
#' @importFrom MAST assay
#' @importFrom stats median
#' @importFrom stats mean
#' @importFrom stats optim
#' 
#' @export
estimation_delta_margin <- function(M) {
	## De-log the gene expression matrix                                     ##
	M <- exp(M)
	p <- nrow(M)
	m <- ncol(M)
	##_______________________________________________________________________##
	## Size factor estimation (adapted from Anders et al. 2010)              ##
	##_______________________________________________________________________##
	## Size factors                                                          ##
	## Pseudo count of 1 to avoid null issues                                ##
	## Computation of the geometric mean of the rows                         ##
	#pp <- apply(M, 1, function(x) prod((x+1)**(1/m)))
	sizeFact <- apply(M, 1, function(x) prod((x+1)**(1/m)))
	#sizeFact <- apply(M, 2, function(x) median(x/pp))
	sizeFact <- apply(M, 2, function(x) median(x/sizeFact))
	##_______________________________________________________________________##
	## Delta estimation (adapted from Yu et al. 2013)                        ##
	##_______________________________________________________________________##
	## access to mu_ik = mean with mu[i], i gene                             ##
	mu <- apply(M, 1, function(x) sum(x/(m*sizeFact)))
	##     var = mean + mean^2 * delta                                       ## 
	## The difference with Yu et al. is that we have access only to the      ##
	## normalized gene expression values                                     ##
	## Naive estimation of variance                                          ##
	#v <- sapply(1:p, function(i) sumSquares(M[i, ]/sizeFact-mu[i])/(m-1))
	est_delta <- sapply(1:p, function(i) sumSquares(M[i, ]/sizeFact-mu[i])/(m-1))
	## Unskrunk dispersion parameter                                         ##
	est_delta <- sapply(1:p,
		#function(i) max(0, (m*v[i]-mu[i]
		function(i) max(0, (m*est_delta[i]-mu[i]
				*sum(1/sizeFact))/(mu[i]**2*sum(1/sizeFact))))#)
	rm(list=c("mu", "sizeFact"))
	weight <- function(zeta) {
		return((sumSquares(est_delta-mean(est_delta))
			/sumSquares(est_delta-zeta))*(p-2)/(p-1))
	}
	## slope of function zeta -> 1/sum_i(est_delta_i-zeta)^2 > 0               ##
	## derivative respect to zeta                                              ##
	slope <- function(zeta) {
		return(-2*sum(est_delta-zeta)/(sumSquares(est_delta-zeta))**2)
	}
	## see reference                                                           ##
	epsilon <- 5e-2
	## Get an unconstrained minimization optimization problem with log-barrier ##
	##         0 > slope(zeta) > -epsilon and zeta > 0                         ##
	## and a penalty on the complexity of the model (L2 norm)                  ##
	f <- function(zeta) {
		return(-log(epsilon+slope(zeta))-log(-slope(zeta)*zeta)+(zeta)**2)
	}
	## Gradient of slope respect to zeta (minus/modulo constants)              ##
	gslope <- function(zeta) {
		return((2*p*sumSquares(est_delta-zeta)
                     -8*(sum(est_delta-zeta)**2))/(sumSquares(est_delta-zeta)**3))
	}
	## Gradient of objective (respect to zeta)                                 ##
	g <- function(zeta) {
		return(-gslope(zeta)/(epsilon+slope(zeta))
			-(slope(zeta)+zeta*gslope(zeta))/(slope(zeta)*zeta)+2*zeta)
	}
	#fest_delta <- sapply(est_delta, f)
	#zeta <- est_delta[which.min(fest_delta)]
	zeta <- est_delta[which.min(sapply(est_delta, f))]
	fit <- optim(par=c(zeta), fn=f, method="BFGS", gr=g)
	zeta <- fit$par
	## Shrinked dispersion parameter (in LOG)                                  ##
	delta <- log((1-weight(zeta))*(est_delta-zeta)+zeta)
	return(delta)
}

#' Estimate of dispersion parameter sigma for copula
#'
#' Estimates of dispersion parameter sigma for copula.
#'
#' @param M data matrix on which estimate sigma
#' @return cvm value of sigma parameter
#' 
#' @importFrom corpcor cor.shrink
#' 
#' @export
estimation_sigma_copula <- function(M) {
	capture.output(cvm <- cor.shrink(t(M))[, ])
	return(cvm)
}

#' Estimate of mixture rates for margins
#'
#' Estimates of mixture rates for margins.
#'
#' @param M data matrix on which estimate mixture rate
#' @param model_dropout model from binary logistic regression
#' @return alpha vector of float values
#' 
#' @importFrom stats predict
#' @importFrom stats mean
#' 
#' @export
estimation_mixture_margin <- function(M, model_dropout) {
	alpha <- sapply(1:nrow(M),
		 function(i) predict(model_dropout, data.frame(Value=median(M[i,])), type="response"))#)
	return(alpha)
}

#' Estimate of mixture rates for copula
#'
#' Estimates of mixture rates for copula.
#'
#' @param M data matrix on which estimate mixture rate
#' @param model_mild models from binary logistic regression
#' @return mix vector of float values
#' 
#' @importFrom stats predict
#' @importFrom stats mean
#' 
#' @export
estimation_mixture_copula <- function(M, replicates, model_mild) {
	replicates <- sapply(replicates, nrow=1, as.numeric)
	mix <- sapply(1:ncol(M), 
		function(j) predict(model_mild[[replicates[j]]], data.frame(Value=mean(M[,j])), type="response"))
	return(mix)
}

######################################
# Mean parameters estimation         #
######################################

#' Estimate mean parameter values
#'
#' Estimates mean parameter values for initialization
#'
#' @param M data matrix on which values are fitted
#' @param peaks peak list associated with @M
#' @param Pj cell pattern for condition #j on which values are fitted
#' @return meanParam estimated margin and copula parameter values
#'
#' @export
estimation_means <- function(M, peaks, Pj, replicates) {
	M <- exp(M)
	m <- ncol(M)
	p <- nrow(M)
	pp <- apply(M, 1, function(x) prod((x+1)**(1/m)))
	sizeFact <- apply(M, 2, function(x) median(x/pp))
	## Anders et al., 2013                               ##
	## Back to log                                       ##
	meanMargin <- log(sapply(1:p, function(i) sum(M[i, ]/(m*sizeFact))))
	## Second highest peak                               ##
	mean1 <- sapply(peaks, function(p) p[2])
	mean1[sapply(mean1, is.na)] <- 0
	## Highest peak                                      ##
	mean2 <- sapply(peaks, function(p) p[1])
	mean1[sapply(mean2, is.na)] <- 0
	## Zero is the baseline gene expression value        ##
	## after the threshold computation by MAST           ##
	selectPeaks <- function(peaks) return(ifelse(length(peaks)==0, 0, mean(peaks)))
	mean1 <- selectPeaks(mean1[Pj==0])
	mean2 <- selectPeaks(mean2[Pj==1])
	## Peaks are already in log                          ##
	meanParam <- c(mean1, mean2, meanMargin)
	return(meanParam)
}

###############################################
## DISTRIBUTION FUNCTIONS                    ##
###############################################

#' Compute cumulative probability for margin
#'
#' Computes cumulative probability for margin which is a Gamma mixture
#' (continuous approximation of Poisson/Negative Binomial mixture). 
#'
#' @param q value vector
#' @param params R named list of parameters
#' @param i index of the margin function in parameter list
#' @param log.p logical for returning log-probability
#' @param lower.tail logical for returning 1-Prob(x <= q) instead of Prob(x <= q)
#' @return p = Prob(x <= q) (unless lower.tail is set to True)
#' 
#' @importFrom stats pgamma
#' @importFrom stats log
#' 
#' @export
pmargin <- function(q, params, i, log.p=FALSE, lower.tail=TRUE) {
	## Same margin parameters for all copulas ##
	margin <- params[[1]]$paramsMargins[[i]]
	y <- margin$alpha*pgamma(q, shape=margin$lambda, scale=1)
		+(1-margin$alpha)*pgamma(q, shape=margin$mu/(1+margin$delta*margin$mu), 
			scale=1+margin$mu*margin$delta)
	if (lower.tail) return(ifelse(log.p, log(y), y))
	else return(ifelse(log.p, log(1-y), 1-y))
}

#' Transform from random variable to uniform random variable
#'
#' Transforms from random variable to uniform random variable
#' ("Probability Integral transform")
#'
#' @param x value vector
#' @param params R named list of parameters
#' @return q transformed vector
#' 
#' @importFrom stats abs
#' 
#' @export
transformQ <- function(x, params) {
	q <- qnorm(as.vector(apply(as.matrix(1:length(x), ncol=1), 1, function(i) {
		pmargin(abs(x[i]), params, i, lower.tail=!(x[i]<0))
	})))
	return(q)
}

#' Compute cumulative probability for copula
#'
#' Computes cumulative probability for copula which is a mixture of
#' Gaussian copulas with Gamma mixture margins
#'
#' @param q transformed value vector
#' @param params R named list of parameters
#' @param j index of the bimodal copula in parameter list
#' @return retvals Prob(C <= q)
#' 
#' @importFrom stats sum
#' 
#' @export
copula <- function(q, params, j) {
	param <- params[[j]]
	retvals <- sum(sapply(1:length(param$mixCop), function(k) {
		param$mixCop[[k]]*pmnormCustom(x=as.numeric(q), mean=param$muCop[[k]], 
		varcov=param$rhoCop[[k]], utrho=param$utrhoCop[[k]])
	}))
	return(retvals)
}

###############################################
## ONE LINER FOR PARAMETER ESTIMATION        ##
###############################################

parameter_estimation <- function(sca, model_dropout, model_mild, peaks, P, options) {
	## Dispersion parameter for each margin/feature ##
	delta <- estimation_delta_margin(assay(sca))
	clock("Estimation of delta")
	## Copula dispersion parameter                  ##
	sigma <- estimation_sigma_copula(assay(sca))
	clock("Estimation of sigma")
	sigma <- nearPDCustom(sigma)
	utrho <- upperTriCustom(sigma)
	clock("Positive definite sigma")
	## Mixture rate parameters                      ##
	alpha <- estimation_mixture_margin(assay(sca), model_dropout)
	clock("Estimation of alpha")
	mix <- estimation_mixture_copula(assay(sca), options$replicates, model_mild)
	clock("Estimation of mix")
	## Copula & Margin mean parameters              ##
	means <- lapply(1:ncol(sca), function(j) estimation_means(assay(sca), peaks, P[, j], options$replicates))
	paramsMargins <- lapply(1:nrow(sca), function(i) list(alpha=alpha[i], lambda=lambda_0, mu=means[[1]][i+2], delta=delta[i]))
	params <- lapply(1:ncol(sca), function(j) {
		muCop <- list(rep(means[[j]][1], nrow(sca)), rep(means[[j]][2], nrow(sca)))
		mixCop <- list(mix[j], 1-mix[j])
		rhoCop <- rep(list(sigma), 2)
		utrhoCop <- rep(list(utrho), 2)
		list(muCop=muCop, mixCop=mixCop, rhoCop=rhoCop, utrhoCop=utrhoCop, paramsMargins=paramsMargins)
	})
	return(params)
}
