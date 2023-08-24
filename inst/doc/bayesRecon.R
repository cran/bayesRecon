## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval=TRUE ### !!!! set to FALSE here to render only the text !!!!
)
set.seed(42)

## ----klippy, echo=FALSE, include=TRUE, eval=FALSE-----------------------------
#  klippy::klippy(position = c('top', 'right'), tooltip_message = 'Copy', tooltip_success = 'Done', color="black")

## ----install, eval=FALSE------------------------------------------------------
#  install.packages('bayesRecon', dependencies = TRUE)

## ----load---------------------------------------------------------------------
library(bayesRecon)

## ----carpart-plot, dpi=300, out.width = "100%", fig.align='center', fig.cap="**Figure 1**: Carpart - monthly car part sales.", fig.dim = c(6, 3)----
layout(mat = matrix(c(1, 2), nrow = 1, ncol = 2), widths = c(2, 1))
plot(carpart, xlab = "Time", ylab = "Car part sales", main = NULL)
hist(carpart, xlab = "Car part sales", main = NULL)

## ----train-test---------------------------------------------------------------
train <- window(carpart, end = c(2001, 3))
test <- window(carpart, start = c(2001, 4))

## ----temp-agg-----------------------------------------------------------------
train.agg <- bayesRecon::temporal_aggregation(train, agg_levels = c(2, 3, 4, 6, 12))
levels <- c("Annual", "Biannual", "4-Monthly", "Quarterly", "2-Monthly", "Monthly")
names(train.agg) <- levels

## ----temp-agg-plot, dpi=300, fig.show="hold", out.width="100%", out.heigth="100%", fig.align='center', fig.cap="**Figure 2**: Carpart - visualization of the aggregated time series.", fig.dim=c(6,3.5)----
par(mfrow = c(2, 3), mai = c(0.6, 0.6, 0.5, 0.5))
for (l in levels) {
  plot(train.agg[[l]], xlab = "Time", ylab = "Car part sales", main = l)
}

## ----hier-fore, cache=TRUE----------------------------------------------------
# install.packages("glarma", dependencies = TRUE)
#library(glarma)

fc.samples <- list()
D <- 20000
fc.count <- 1
# iterating over the temporal aggregation levels
for (l in seq_along(train.agg)) {
  f.level <- frequency(train.agg[[l]])
  print(paste("Forecasting at ", levels[l], "...", sep = ""))
  # fit an independent model for each aggregation level
  model <- glarma::glarma(
    train.agg[[l]],
    phiLags = if (f.level == 1) 1 else 1:(min(6, f.level - 1)),
    thetaLags = if (f.level == 1) NULL else f.level,
    X = cbind(intercept = rep(1, length(train.agg[[l]]))),
    offset = cbind(intercept = rep(0, length(train.agg[[l]]))),
    type = "Poi"
  )
  # forecast 1 year ahead
  h <- f.level
  tmp <- matrix(data = NA, nrow = h, ncol = D)
  for (s in (1:D)) {
    # each call to 'forecast.glarma' returns a simulation path
    tmp[, s] <- glarma::forecast(
      model,
      n.ahead = h,
      newdata = cbind(intercept = rep(1, h)),
      newoffset = rep(0, h)
    )$Y
  }
  # collect the forecasted samples
  for (i in 1:h) {
    fc.samples[[fc.count]] <- tmp[i, ]
    fc.count <- fc.count + 1
  }
}


## ----S------------------------------------------------------------------------
recon.matrices <- bayesRecon::get_reconc_matrices(agg_levels = c(2, 3, 4, 6, 12), h = 12)
# Summing matrix
S <- recon.matrices$S 
# A matrix
A <- recon.matrices$A 

## ----reconc-------------------------------------------------------------------
recon.res <- bayesRecon::reconc_BUIS(
  S,
  base_forecasts = fc.samples,
  in_type = "samples",
  distr = "discrete",
  seed = 42
)

## ----res----------------------------------------------------------------------
reconciled_samples <- recon.res$reconciled_samples
dim(reconciled_samples)


## ----metrics------------------------------------------------------------------
# install.packages("scoringRules", dependencies = TRUE)
library(scoringRules)

ae.fc <- list()
ae.reconc <- list()
crps.fc <- list()
crps.reconc <- list()
for (h in 1:length(test)) {
  y.hat_ <- median(fc.samples[[nrow(A) + h]])
  y.reconc_ <- median(recon.res$bottom_reconciled_samples[, h])
  # Compute Absolute Errors
  ae.fc[[h]] <- abs(test[h] - y.hat_)
  ae.reconc[[h]] <- abs(test[h] - y.reconc_)
  # Compute Continuous Ranked Probability Score (CRPS)
  crps.fc[[h]] <-
    scoringRules::crps_sample(y = test[h], dat = fc.samples[[nrow(A) + h]])
  crps.reconc[[h]] <-
    scoringRules::crps_sample(y = test[h], dat = recon.res$bottom_reconciled_samples[, h])
}

mae.fc <- mean(unlist(ae.fc))
mae.reconc <- mean(unlist(ae.reconc))
crps.fc <- mean(unlist(crps.fc))
crps.reconc <- mean(unlist(crps.reconc))
metrics <- data.frame(
  row.names = c("MAE", "CRPS"),
  base.forecasts = c(mae.fc, crps.fc),
  reconciled.forecasts = c(mae.reconc, crps.reconc)
)
metrics


## ----m3-plot, dpi=300, out.width = "100%", fig.align='center', fig.cap="**Figure 3**: M3 - N1485 time series.", fig.dim = c(6, 3)----
plot(M3example$train, xlab = "Time", ylab = "y", main = "N1485")

## ----m3-agg-------------------------------------------------------------------
train.agg <- bayesRecon::temporal_aggregation(M3example$train)
levels <- c("Annual", "Biannual", "4-Monthly", "Quarterly", "2-Monthly", "Monthly")
names(train.agg) <- levels

## ----m3-fore------------------------------------------------------------------
# install.packages("forecast", dependencies = TRUE)
library(forecast)

H <- length(M3example$test)
H

fc <- list()
level.idx <- 1
fc.idx <- 1
for (level in train.agg) {
  level.name <- names(train.agg)[level.idx]
  # fit an ETS model for each temporal level
  model <- ets(level)
  # generate forecasts for each level within 18 months
  h <- floor(H / (12 / frequency(level)))
  print(paste("Forecasting at ", level.name, ", h=", h, "...", sep = ""))
  level.fc <- forecast(model, h = h)
  # save mean and sd of the gaussian predictive distribution
  for (i in 1:h) {
    fc[[fc.idx]] <- c(level.fc$mean[[i]],
                      (level.fc$upper[, "95%"][[i]] - level.fc$mean[[i]]) / qnorm(0.975))
    fc.idx <- fc.idx + 1
  }
  level.idx <- level.idx + 1
}

## ----m3-rmat, dpi=300, out.width = '70%', fig.align='center', fig.cap="**Figure 4**: M3 - The aggregating matrix A (red=1, yellow=0).", fig.dim = c(8, 8)----
rmat <- get_reconc_matrices(agg_levels = c(2, 3, 4, 6, 12), h = 18)

par(mai = c(1,1,0.5,0.5))
image(1:ncol(rmat$A), 1:nrow(rmat$A), 
      t(apply(t(rmat$A),1,rev)), 
      xaxt='n', yaxt='n', ylab = "", xlab = levels[6])
axis(1, at=1:ncol(rmat$A), label=1:ncol(rmat$A), las=1)
axis(2, at=c(23,22,19,15,9), label=levels[1:5], las=2)

## ----m3-reco------------------------------------------------------------------
recon.gauss <- bayesRecon::reconc_gaussian(
  S = rmat$S,
  base_forecasts.mu = sapply(fc, "[[", 1),
  base_forecasts.Sigma = diag(sapply(fc, "[[", 2)) ^ 2
)

reconc.buis <- bayesRecon::reconc_BUIS(
  S = rmat$S,
  base_forecasts = fc,
  in_type = "params",
  distr = "gaussian",
  num_samples = 20000,
  seed = 42
)

# check that the algorithms return consistent results
round(rbind(
  c(recon.gauss$upper_reconciled_mean, recon.gauss$bottom_reconciled_mean),
  rowMeans(reconc.buis$reconciled_samples)
))

## ----m3-plotfore, dpi=300, out.width = "100%", fig.align='center', fig.cap="**Figure 5**: M3 - visualization of base and reconciled forecasts. The black line is the actual data (dashed is the test). The orange line is the forecasted mean, the blu line the reconciled mean. Shadow regions show the 95% prediction intervals.", fig.dim = c(6, 4)----
yhat.mu <- tail(sapply(fc, "[[", 1), 18)
yhat.sigma <- tail(sapply(fc, "[[", 2), 18)
yhat.hi95 <- qnorm(0.975, mean = yhat.mu, sd = yhat.sigma)
yhat.lo95 <- qnorm(0.025, mean = yhat.mu, sd = yhat.sigma)
yreconc.mu <- rowMeans(reconc.buis$bottom_reconciled_samples)
yreconc.hi95 <- apply(reconc.buis$bottom_reconciled_samples, 1, 
                      function(x) quantile(x, 0.975))
yreconc.lo95 <- apply(reconc.buis$bottom_reconciled_samples, 1, 
                      function(x) quantile(x, 0.025))

ylim_ <- c(min(M3example$train, M3example$test, yhat.lo95, yreconc.lo95) - 1, 
           max(M3example$train, M3example$test, yhat.hi95, yreconc.hi95) + 1)

plot(M3example$train, xlim = c(1990, 1995.6), ylim = ylim_, 
     ylab = "y", main = "N1485 Forecasts")
lines(M3example$test, lty = "dashed")
lines(ts(yhat.mu, start = start(M3example$test), frequency = 12), 
      col = "coral", lwd = 2)
lines(ts(yreconc.mu, start = start(M3example$test), frequency = 12), 
      col = "blue2", lwd = 2)
xtest <- time(M3example$test)
polygon(c(xtest, rev(xtest)), c(yhat.mu, rev(yhat.hi95)), 
        col = "#FF7F5066", border = "#FF7F5066")
polygon(c(xtest, rev(xtest)), c(yhat.mu, rev(yhat.lo95)), 
        col = "#FF7F5066", border = "#FF7F5066")
polygon(c(xtest, rev(xtest)), c(yreconc.mu, rev(yreconc.hi95)), 
        col = "#0000EE4D", border = "#0000EE4D")
polygon(c(xtest, rev(xtest)), c(yreconc.mu, rev(yreconc.lo95)), 
        col = "#0000EE4D", border = "#0000EE4D")

## ----infants-forecasts--------------------------------------------------------
# install.packages("forecast", dependencies = TRUE)
library(forecast)

fc <- list()
residuals <- matrix(NA,
                    nrow = length(infantMortality$total),
                    ncol = length(infantMortality))
fc.idx <- 1
for (s in infantMortality) {
  s.name <- names(infantMortality)[fc.idx]
  print(paste("Forecasting at ", s.name, "...", sep = ""))
  # fit an auto.arima model and forecast with h=1
  model <- auto.arima(s)
  s.fc <- forecast(model, h = 1)
  # save mean and sd of the gaussian predictive distribution
  fc[[s.name]] <- c(s.fc$mean,
                    (s.fc$upper[, "95%"][[1]] - s.fc$mean) / qnorm(0.975))
  residuals[, fc.idx] <- s.fc$residuals
  fc.idx <- fc.idx + 1
}

## ----infants-s, dpi=300, out.width = '70%', fig.align='center', fig.cap="**Figure 6**: Infants mortality - The aggregating matrix A (red=1, yellow=0).", fig.dim = c(8, 8)----
# we have 16 bottom time series, and 11 upper time series
A <- matrix(data = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                     0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,
                     0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,
                     0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,
                     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,
                     1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,
                     0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1), byrow=TRUE, ncol = 16)
S <- rbind(A, diag(16))

# plot of A
par(mai = c(1.5,1,0.5,0.5))
image(1:ncol(A), 1:nrow(A), 
      t(apply(t(A),1,rev)), 
      xaxt='n', yaxt='n', ann=FALSE)
axis(1, at=1:ncol(A), label=names(infantMortality)[12:27], las=2)
axis(2, at=c(1:11), label=rev(names(infantMortality)[1:11]), las=2)

## ----infants reconc-----------------------------------------------------------
# means
mu <- sapply(fc, "[[", 1)
# Shrinkage covariance
shrink.res <- bayesRecon::schaferStrimmer_cov(residuals)
print(paste("The estimated shrinkage intensity is", round(shrink.res$lambda_star, 3)))
Sigma <- shrink.res$shrink_cov

## ----infants-recon------------------------------------------------------------
recon.gauss <- bayesRecon::reconc_gaussian(S,
                                           base_forecasts.mu = mu,
                                           base_forecasts.Sigma = Sigma)

# check coherence
(A %*% recon.gauss$bottom_reconciled_mean) - recon.gauss$upper_reconciled_mean

