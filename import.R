################################################################
#                                                              #
#               MineClusters OPTIONS & LIBRARIES & FILES       #
#                                                              #
################################################################

#####################
# Libraries         #
#####################

## Import functions for multinormal distributions       ##
library(mnormt)
## Import function to create folds for cross-validation ##
library(caret)
## Import function to compute gene expression threshold ##
library(MAST)
## Parallelization library                              ##
library(foreach)
library(doParallel)
library(doRNG)
## Hierarchical clustering                              ##
library(fastcluster)
## Import function to estimate covariance matrix        ##
library(corpcor)
library(stats)
library(bigstatsr)

#####################
# Options           #
#####################

options(expressions = 5e5, warn=-1)
nacores = parallel::detectCores()
ncores <- ifelse(is.null(nacores), 1, ifelse(nacores<3, 1, min(nacores-2, 5)))
lambda_0 = 0.1

#####################
# Files             #
#####################

source("utils/utils.R")
source("utils/pattern.R")
source("utils/computation.R")
source("utils/feature_selection.R")
