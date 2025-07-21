## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(bayesRecon)

## ----M5hier, fig.cap="**Figure 1**: graph of the M5 hierarchy.", out.width = '100%', echo = FALSE----
knitr::include_graphics("img/M5store_hier.png")

## ----InitializeHierarchy------------------------------------------------------
# Hierarchy composed by 3060 time series: 3049 bottom and 11 upper
n_b <- 3049
n_u <- 11
n <- n_b + n_u   

# Load matrix A 
A <- M5_CA1_basefc$A

# Load base forecasts:
base_fc_upper  <- M5_CA1_basefc$upper
base_fc_bottom <- M5_CA1_basefc$bottom

# We will save all the results in the list rec_fc
rec_fc <- list(
            Gauss      = list(),
            Mixed_cond = list(),
            TD_cond    = list()
          )

## ----readBaseForecasts--------------------------------------------------------
# Parameters of the upper base forecast distributions
mu_u <- unlist(lapply(base_fc_upper, "[[", "mu"))  # upper means
# Compute the (shrinked) covariance matrix of the residuals
residuals.upper <- lapply(base_fc_upper, "[[", "residuals")
residuals.upper <- t(do.call("rbind", residuals.upper))
Sigma_u <- schaferStrimmer_cov(residuals.upper)$shrink_cov
  
# Parameters of the bottom base forecast distributions
mu_b <- c()
sd_b <- c()
for (fc_b in base_fc_bottom) {
  pmf <- fc_b$pmf
  mu_b <- c(mu_b, PMF.get_mean(pmf))
  sd_b <- c(sd_b, PMF.get_var(pmf)**0.5)
}
Sigma_b <- diag(sd_b**2)

# Mean and covariance matrix of the base forecasts
base_forecasts.mu <- c(mu_u,mu_b)
base_forecasts.Sigma <- matrix(0, nrow = n, ncol = n)
base_forecasts.Sigma[1:n_u,1:n_u] <- Sigma_u
base_forecasts.Sigma[(n_u+1):n,(n_u+1):n] <- Sigma_b

## ----GaussianReconciliation---------------------------------------------------
# Gaussian reconciliation 
start <- Sys.time()       
gauss <- reconc_gaussian(A, base_forecasts.mu, base_forecasts.Sigma)
stop <- Sys.time()

rec_fc$Gauss <- list(mu_b    = gauss$bottom_reconciled_mean,
                     Sigma_b = gauss$bottom_reconciled_covariance,
                     mu_u    = A %*% gauss$bottom_reconciled_mean,
                     Sigma_u = A %*% gauss$bottom_reconciled_covariance %*% t(A))

Gauss_time <- as.double(round(difftime(stop, start, units = "secs"), 2))
cat("Time taken by Gaussian reconciliation: ", Gauss_time, "s")

## ----MixedCondReconciliation--------------------------------------------------
seed <- 1
N_samples_IS <- 5e4

# Base forecasts
fc_upper_4rec <- list(mu=mu_u, Sigma=Sigma_u)
fc_bottom_4rec <- lapply(base_fc_bottom, "[[", "pmf")  # list of PMFs

# MixCond reconciliation
start <- Sys.time()       
mix_cond <- reconc_MixCond(A, fc_bottom_4rec, fc_upper_4rec, bottom_in_type = "pmf",
                          num_samples = N_samples_IS, return_type = "pmf", seed = seed)
stop <- Sys.time()       

rec_fc$Mixed_cond <- list(
  bottom = mix_cond$bottom_reconciled$pmf,
  upper  = mix_cond$upper_reconciled$pmf,
  ESS    = mix_cond$ESS
  )

MixCond_time <- as.double(round(difftime(stop, start, units = "secs"), 2))
cat("Computational time for Mix-cond reconciliation: ", MixCond_time, "s")

## ----TDcondReconciliation-----------------------------------------------------
N_samples_TD <- 1e4

# TDcond reconciliation
start <- Sys.time()     
td <- reconc_TDcond(A, fc_bottom_4rec, fc_upper_4rec,
                   bottom_in_type = "pmf", num_samples = N_samples_TD, 
                   return_type = "pmf", seed = seed)
stop <- Sys.time()  

## ----print computational time TD cond-----------------------------------------
rec_fc$TD_cond <- list(
  bottom = td$bottom_reconciled$pmf,
  upper  = td$upper_reconciled$pmf
  )

TDCond_time <- as.double(round(difftime(stop, start, units = "secs"), 2))
cat("Computational time for TD-cond reconciliation: ", TDCond_time, "s")

## ----InitializeMetrics--------------------------------------------------------
# Parameters for computing the scores
alpha <- 0.1   # MIS uses 90% coverage intervals
jitt <- 1e-9   # jitter for numerical stability 

# Save actual values
actuals_u <- unlist(lapply(base_fc_upper, "[[", "actual"))
actuals_b <- unlist(lapply(base_fc_bottom, "[[", "actual"))
actuals <- c(actuals_u, actuals_b)

# Scaling factor for computing MASE
Q_u <- M5_CA1_basefc$Q_u
Q_b <- M5_CA1_basefc$Q_b
Q   <- c(Q_u, Q_b)

# Initialize lists to save the results
mase <- list()
mis  <- list()
rps  <- list()

## ----metricFunctions,include=FALSE--------------------------------------------
# Functions for computing the scores of a PMF
AE_pmf <- function(pmf, actual) {
  return(abs(PMF.get_quantile(pmf,p=0.5) - actual))
}
MIS_pmf <- function(pmf, actual, alpha) {
  u <- PMF.get_quantile(pmf, p=1-(alpha/2))
  l <- PMF.get_quantile(pmf, p=alpha/2)
  return(u - l + (2/alpha)*(l - actual)*(actual < l) + 
                 (2/alpha)*(actual - u)*(actual > u) )
}
RPS_pmf <- function(pmf, actual) {
  cdf <- cumsum(pmf) / sum(pmf)
  M <- length(cdf)
  # if actual is outside the supp of pmf, add ones to the end of the cdf:
  if (actual >= M) {  
    cdf <- c(cdf, rep(1, (actual-M+1)))
    M <- length(cdf)
  }
  cdf_act <- (0:(M-1)) >= actual   # unit step function in actual
  crps_ <- sum((cdf - cdf_act)**2)
  return(crps_)
}

# Function for computing the MIS of (possibly truncated) Gaussian forecasts
MIS_gauss <- function(mus, sds, actuals, alpha, trunc=FALSE) {
  z <- qnorm(1-(alpha/2))
  u <- mus + z * sds
  l <- mus - z * sds
  # If it is truncated, set negative quantiles to zero
  if (trunc) {
    l[l<0] <- 0
    u[u<0] <- 0
  }  
  return(u - l + (2/alpha)*(l - actuals)*(actuals < l) + 
                 (2/alpha)*(actuals - u)*(actuals > u) )
}

## ----computeScores------------------------------------------------------------
# Compute scores for the base forecasts
# Upper
mu_u <- unlist(lapply(base_fc_upper, "[[", "mu"))
sd_u <- unlist(lapply(base_fc_upper, "[[", "sigma"))
mase$base[1:n_u] <- abs(mu_u - actuals_u) / Q_u
mis$base[1:n_u]  <- MIS_gauss(mu_u, sd_u, actuals_u, alpha)
rps$base[1:n_u]  <- scoringRules::crps(actuals_u, "norm", mean=mu_u, sd=sd_u)
# Bottom
pmfs = lapply(base_fc_bottom, "[[", "pmf")
mase$base[(n_u+1):n] <- mapply(AE_pmf, pmfs, actuals_b) / Q_b
mis$base[(n_u+1):n]  <- mapply(MIS_pmf, pmfs, actuals_b, MoreArgs = list(alpha=alpha))
rps$base[(n_u+1):n]  <- mapply(RPS_pmf, pmfs, actuals_b)

# Compute scores for Gauss reconciliation
mu <- c(rec_fc$Gauss$mu_u, rec_fc$Gauss$mu_b)
sd <- c(diag(rec_fc$Gauss$Sigma_u), diag(rec_fc$Gauss$Sigma_b))**0.5
sd <- sd + jitt
mase$Gauss <- abs(mu - actuals) / Q
mis$Gauss <- MIS_gauss(mu, sd, actuals, alpha)
rps$Gauss <- scoringRules::crps(actuals, "norm", mean=mu, sd=sd)

# Compute scores for Mix-cond reconciliation
pmfs <- c(rec_fc$Mixed_cond$upper, rec_fc$Mixed_cond$bottom)
mase$MixCond <- mapply(AE_pmf, pmfs, actuals) / Q
mis$MixCond  <- mapply(MIS_pmf, pmfs, actuals, MoreArgs = list(alpha=alpha))
rps$MixCond  <- mapply(RPS_pmf, pmfs, actuals)

# Compute scores for TD-cond reconciliation
pmfs <- c(rec_fc$TD_cond$upper, rec_fc$TD_cond$bottom)
mase$TDcond <- mapply(AE_pmf, pmfs, actuals) / Q
mis$TDcond  <- mapply(MIS_pmf, pmfs, actuals, MoreArgs = list(alpha=alpha))
rps$TDcond  <- mapply(RPS_pmf, pmfs, actuals)

## ----skillScoreFunction,include=FALSE-----------------------------------------
# Function for computing the skill score
skill.score <- function(ref, met) {
  s <- (2 * (ref - met) / (ref + met)) * 100
  s[is.na(s)] <- 0  # if both numerator and denominator are 0, set skill score to 0
  return(s)
}

## ----ComputeScores------------------------------------------------------------
scores <- list(
        mase = mase,
        mis  = mis,
        rps = rps
      )
scores_ = names(scores)

ref_met  <- "base"
methods_ <- c("Gauss", "MixCond", "TDcond")

# For each score and method we compute the skill score with respect to the base forecasts
skill_scores <- list()
for (s in scores_) {
  skill_scores[[s]] <- list()
  for (met in methods_) {
    skill_scores[[s]][["upper"]][[met]] <- skill.score(scores[[s]][[ref_met]][1:n_u], 
                                                       scores[[s]][[met]][1:n_u])
    skill_scores[[s]][["bottom"]][[met]] <- skill.score(scores[[s]][[ref_met]][(n_u+1):n], 
                                                        scores[[s]][[met]][(n_u+1):n])
    }
}

## ----AverageScores------------------------------------------------------------
mean_skill_scores <- list()

for (s in scores_) {
  mean_skill_scores[[s]] <- rbind(data.frame(lapply(skill_scores[[s]][["upper"]], mean)),
                                    data.frame(lapply(skill_scores[[s]][["bottom"]], mean))
                                     )
  rownames(mean_skill_scores[[s]]) <- c("upper","bottom")
}

## ----printMASETable-----------------------------------------------------------
knitr::kable(mean_skill_scores$mase,digits = 2,caption = "Mean skill score on MASE.",align = 'lccc')

## ----PrintMIStable------------------------------------------------------------
knitr::kable(mean_skill_scores$mis,digits = 2,caption = "Mean skill score on MIS.")

## ----printRPStable------------------------------------------------------------
knitr::kable(mean_skill_scores$rps,digits = 2,caption = "Mean skill score on RPS.")

## ----MASEboxplots, fig.cap="**Figure 2**: boxplot of MASE skill scores for upper and bottom time series.", fig.width=7,fig.height=8----
custom_colors <- c("#a8a8e4", 
                  "#a9c7e4",
                  "#aae4df")

# Boxplots of MASE skill scores
par(mfrow = c(2, 1))
boxplot(skill_scores$mase$upper, main = "MASE upper time series", 
        col = custom_colors, ylim = c(-80,80))
abline(h=0,lty=3)
boxplot(skill_scores$mase$bottom, main = "MASE bottom time series", 
        col = custom_colors, ylim = c(-200,200))
abline(h=0,lty=3)

## ----setupParams,eval=TRUE,include=FALSE--------------------------------------
par(mfrow = c(1, 1))

## ----MISboxplots, fig.cap="**Figure 3**: boxplot of MIS skill scores for upper and bottom time series.", fig.width=7,fig.height=8----
# Boxplots of MIS skill scores
par(mfrow = c(2, 1))
boxplot(skill_scores$mis$upper, main = "MIS upper time series", 
        col = custom_colors, ylim = c(-150,150))
abline(h=0,lty=3)
boxplot(skill_scores$mis$bottom, main = "MIS bottom time series", 
        col = custom_colors, ylim = c(-200,200))
abline(h=0,lty=3)

## ----eval=TRUE,include=FALSE--------------------------------------------------
par(mfrow = c(1, 1))

## ----RPSboxplots,fig.cap="**Figure 4**: boxplot of RPS skill scores for upper and bottom time series.", fig.width=7,fig.height=8----
# Boxplots of RPS skill scores
par(mfrow = c(2,1))
boxplot(skill_scores$rps$upper, main = "RPS upper time series", 
        col = custom_colors, ylim = c(-80,80))
abline(h=0,lty=3)
boxplot(skill_scores$rps$bottom, main = "RPS bottom time series", 
        col = custom_colors, ylim = c(-200,200))
abline(h=0,lty=3)

## ----eval=TRUE,include=FALSE--------------------------------------------------
par(mfrow = c(1, 1))

