#' causalGBN: a package for estimating causal Gaussian Bayesian Networks.
#'
#' The packages allows to estimate causal Gaussian Bayesian Networks
#'  from a mixture of observation and intervention experiment.
#'  When the DAG structure G is known, parameters are estimated through
#'  simple linear regression.
#'  In the case where the DAG G is unknown, the posterior distribution
#'  P(G | data) is derived through a Metropolis-Hasting algorithm (MC3)
#'  where the score of each DAG G is obtained through a
#'  BIC (eBIC by default).
#'
#' @section Main functions:
#' The main functions ...
#'
#' @docType package
#' @name causalGBN
NULL
