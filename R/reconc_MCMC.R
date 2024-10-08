#-------------------------------------------------------------------------------
#' @title MCMC for Probabilistic Reconciliation of forecasts via conditioning
#'
#' @description
#'
#' Uses Markov Chain Monte Carlo algorithm to draw samples from the reconciled
#' forecast distribution, which is obtained via conditioning.
#'
#' This is a bare-bones implementation of the Metropolis-Hastings algorithm, 
#' we suggest the usage of tools to check the convergence.
#' The function only works with Poisson or Negative Binomial base forecasts.
#'
#' The function [reconc_BUIS()] is generally faster on most hierarchies.
#'
#'
#'
#' @param A aggregation matrix (n_upper x n_bottom).
#' @param base_forecasts list of the parameters of the base forecast distributions, see details.
#' @param distr a string describing the type of predictive distribution.
#' @param num_samples number of samples to draw using MCMC.
#' @param tuning_int number of iterations between scale updates of the proposal.
#' @param init_scale initial scale of the proposal.
#' @param burn_in number of initial samples to be discarded.
#' @param seed seed for reproducibility.
#'
#' @details
#'
#' The parameter `base_forecast` is a list containing n = n_upper + n_bottom elements.
#' Each element is a list containing the estimated:
#'
#' * mean and sd for the Gaussian base forecast, see \link[stats]{Normal}, if `distr`='gaussian';
#' * lambda for the Poisson base forecast, see \link[stats]{Poisson}, if `distr`='poisson';
#' * size and prob (or mu) for the negative binomial base forecast, see \link[stats]{NegBinomial}, if `distr`='nbinom'.
#'
#' The first n_upper elements of the list are the upper base forecasts, in the order given by the rows of A.
#' The elements from n_upper+1 until the end of the list are the bottom base forecasts, in the order given by the columns of A.
#'
#' @return A list containing the reconciled forecasts. The list has the following named elements:
#'
#' * `bottom_reconciled_samples`: a matrix (n_bottom x `num_samples`) containing reconciled samples for the bottom time series;
#' * `upper_reconciled_samples`: a matrix (n_upper x `num_samples`) containing reconciled samples for the upper time series;
#' * `reconciled_samples`: a matrix (n x `num_samples`) containing the reconciled samples for all time series.
#'
#' @examples
#'
#'library(bayesRecon)
#'
#'# Create a minimal hierarchy with 2 bottom and 1 upper variable
#'rec_mat <- get_reconc_matrices(agg_levels=c(1,2), h=2)
#'A <- rec_mat$A
#'
#'#Set the parameters of the Poisson base forecast distributions
#'lambda1 <- 2
#'lambda2 <- 4
#'lambdaY <- 9
#'lambdas <- c(lambdaY,lambda1,lambda2)
#'
#'base_forecasts = list()
#'for (i in 1:length(lambdas)) {
#'  base_forecasts[[i]] = list(lambda = lambdas[i])
#'}
#'
#'#Sample from the reconciled forecast distribution using MCMC
#'mcmc = reconc_MCMC(A, base_forecasts, distr = "poisson",
#'                   num_samples = 30000, seed = 42)
#'samples_mcmc <- mcmc$reconciled_samples
#'
#'#Compare the reconciled means with those obtained via BUIS
#'buis = reconc_BUIS(A, base_forecasts, in_type="params",
#'                    distr="poisson", num_samples=100000, seed=42)
#'samples_buis <- buis$reconciled_samples
#'
#'print(rowMeans(samples_mcmc))
#'print(rowMeans(samples_buis))
#'
#' @references
#' Corani, G., Azzimonti, D., Rubattu, N. (2024). 
#' *Probabilistic reconciliation of count time series*. 
#' International Journal of Forecasting 40 (2), 457-469.
#' \doi{10.1016/j.ijforecast.2023.04.003}.
#'
#' @seealso
#' [reconc_BUIS()]
#'
#' @export
reconc_MCMC <- function(A,
                        base_forecasts,
                        distr,
                        num_samples = 10000,
                        tuning_int = 100,
                        init_scale = 1,
                        burn_in = 1000,
                        seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # Ensure that data inputs are valid
  if (distr == "gaussian") {
    stop("MCMC for Gaussian distributions is not implemented")
  }
  
  n_bottom <- ncol(A)
  n_ts <- nrow(A)+ncol(A)
  
  # Transform distr into list
  if (!is.list(distr)) {
    distr = rep(list(distr), n_ts)
  }
  # Check input
  .check_input_BUIS(A, base_forecasts, in_type = as.list(rep("params", n_ts)), distr = distr)

  # the first burn_in samples will be removed
  num_samples <- num_samples + burn_in

  # Set the covariance matrix of the proposal (identity matrix)
  cov_mat_prop <- diag(rep(1,n_bottom))

  #Set the counter for tuning
  c_tuning <- tuning_int

  # Initialize acceptance counter
  accept_count <- 0

  # Create empty matrix for the samples from MCMC
  b <- matrix(nrow = num_samples, ncol = n_bottom)

  # Get matrix A and bottom base forecasts
  split_hierarchy.res <- list(
    A = A,
    upper = base_forecasts[1:nrow(A)],
    bottom = base_forecasts[(nrow(A)+1):n_ts],
    upper_idxs = 1:nrow(A),
    bottom_idxs = (nrow(A)+1):n_ts
  )
  bottom_base_forecasts <- split_hierarchy.res$bottom
  bottom_distr <- distr[split_hierarchy.res$bottom_idxs]

  # Initialize first sample (draw from base distribution)
  b[1,] <- .initialize_b(bottom_base_forecasts, bottom_distr)

  # Initialize prop list
  old_prop <- list(
    "b" = b[1,],
    "scale" = init_scale
  )

  # Run the chain
  for (i in 2:num_samples) {

    if (c_tuning == 0) {
      old_prop$acc_rate <- accept_count / tuning_int  #set acc_rate
      accept_count <- 0                               #reset acceptance counter
      c_tuning <- tuning_int                          #reset tuning counter
    }

    prop <- .proposal(old_prop, cov_mat_prop)
    b_prop <- prop$b
    alpha <- .accept_prob(b_prop, b[i-1,], A, distr, base_forecasts)

    if (stats::runif(1) < alpha) {
      b[i,] <- b_prop
      accept_count <- accept_count + 1
    } else {
      b[i,] <- b[i-1,]
    }

    old_prop <- list(
      "b" = b[i,],
      "scale" = prop$scale
    )

    c_tuning <- c_tuning - 1

  }

  b_samples <- t(b[(burn_in+1) : num_samples, ]) #output shape: n_bottom x num_samples
  u_samples <- A %*% b_samples
  y_samples <- rbind(u_samples,b_samples)

  out = list(
    bottom_reconciled_samples = b_samples,
    upper_reconciled_samples = u_samples,
    reconciled_samples = y_samples
  )

  return(out)

}


#-------------------------------------------------------------------------------
##################################
.initialize_b <- function(bottom_base_forecasts, bottom_distr) {

  b <- c()
  for (i in 1:length(bottom_distr)) {
    b[i] <- .distr_sample(bottom_base_forecasts[[i]], bottom_distr[[i]], 1)
  }

  return(b)

}



#-------------------------------------------------------------------------------
##################################
# @title Compute acceptance probability
# @param b proposal state
# @param b0 current state
# @param A aggregation matrix
# @param distr list of strings specifying the distribution of each variable
# @param params list of the parameters of the distributions
# @return the acceptance probability alpha
.accept_prob <- function(b, b0, A, distr, params) {

  alpha <- .target_pmf(b, A, distr, params) / .target_pmf(b0, A, distr, params)

  return(min(1,alpha))

}


##################################
.target_pmf <- function(b, A, distr, params) {

  n_ts <- nrow(A)+ncol(A)

  y <- rbind(A,diag(nrow=ncol(A),ncol=ncol(A))) %*% b

  pmf <- 1
  for (j in 1:n_ts) {
    pmf <- pmf * .distr_pmf(y[[j]], params[[j]], distr[[j]])
  }

  return(pmf)

}


#-------------------------------------------------------------------------------
##################################
# @title Generate a proposal for MCMC step
# @param prev_prop a list containing b, scale, acc_rate
# @param cov_mat_prop the covariance matrix (ncol(cov_mat_prop)=length(b)) for the normal proposal
# @return a list containing the new b, scale, acc_rate
.proposal <- function(prev_prop, cov_mat_prop){

  # extract previous proposal
  b0 <- prev_prop$b
  old_scale <- prev_prop$scale
  acc_rate <- prev_prop$acc_rate

  if(!is.null(acc_rate)){
    # Compute the scaling parameter
    scale <- .tune(old_scale,acc_rate)
  }else{
    scale <- old_scale
  }


  n_x <- length(b0)

  if(n_x != ncol(cov_mat_prop)){
    stop(sprintf("Error in .proposal: previous state dim (%s) and
                 covariance dim (%s) do not match", n_x, ncol(cov_mat_prop)))
  }

  # We always do an independent Gaussian proposal
  dd <- stats::rnorm(n_x)*diag(cov_mat_prop)*scale

  # CHANGE HERE for mixed case
  dd <- round(dd,0)
  b = b0+dd

  prop <- list(b=b, acc_rate=acc_rate, scale=scale)
  return(prop)
}


# scaling parameter
# we use the same variance adaptation used in pymc3
# Rate    Variance adaptation
# ----    -------------------
# <0.001        x 0.1
# <0.05         x 0.5
# <0.2          x 0.9
# >0.5          x 1.1
# >0.75         x 2
# >0.95         x 10
.tune <- function(scale, acc_rate){
  if(acc_rate<0.001){
    return(scale*0.1)
  }else if(acc_rate<0.05){
    return(scale*0.5)
  }else if(acc_rate<0.2){
    return(scale*0.9)
  }else if(acc_rate>0.5){
    return(scale*1.1)
  }else if(acc_rate>0.75){
    return(scale*2)
  }else if(acc_rate>0.95){
    return(scale*10)
  }

  return(scale)

}

