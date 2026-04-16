## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load, message=FALSE, warning=FALSE---------------------------------------
# load the packages
library(bayesRecon)
library(forecast) # base forecasts
library(ggplot2)  # plots

## ----swiss-tourism-plot, dpi=300, out.width = "100%", fig.align='center', fig.cap="**Figure 1**: Swiss tourism - monthly Swiss overnight stays.", fig.dim = c(8, 4)----
# save all time series 
Y = swiss_tourism$ts

# plot the first (top) time series
autoplot(Y[,1], ylab = "Overnight stays in Switzerland",linewidth=0.9)+
  scale_y_continuous(labels = function(x) paste0(formatC(x / 1e6, format = "g"), "M"))

## ----Swiss tourism------------------------------------------------------------
# Save aggregation matrix
A = swiss_tourism$agg_mat

# Number of bottom and upper time series
n_b = ncol(A)
n_u = nrow(A)
n = n_b + n_u

# Frequency is monthly:
print(frequency(Y))

# Select the length of the training set and the forecast horizon
L = 60    
h = 1    

# Select the training set and the actuals for the forecast horizon
train = window(Y, end = time(Y)[L])
actuals = window(Y, start= time(Y)[L + 1], end = time(Y)[L + h])

## ----Base-forecasts and residuals computation---------------------------------
# Compute base forecasts and residuals for each time series
base_fc = rep(NA, n)
res = matrix(NA, ncol = n, nrow = L)
for (i in 1:n){
    models = forecast::ets(train[,i], model = "AZZ")
    f = forecast(models, h = h)
    base_fc[i] = f$mean
    res[, i] = models$residuals
}

## ----t-reconciliation---------------------------------------------------------
t_rec_results = reconc_t(A, base_fc_mean = base_fc, 
                         y_train = train, 
                         residuals = res,
                         return_parameters = TRUE,
                         return_upper = TRUE)

## ----Base---------------------------------------------------------------------
# Base forecasts
Base_mean = base_fc
Base_cov_mat = crossprod(res)/nrow(res) # covariance of the residuals

## ----Gaussian-----------------------------------------------------------------
# Gaussian/MinT: compute reconciliation with bayesRecon
gauss_results = reconc_gaussian(A, base_fc, residuals = res, return_upper = TRUE)
# Reconciled mean for the whole hierarchy:
MinT_reconciled_mean = c(gauss_results$upper_rec_mean, 
                         gauss_results$bottom_rec_mean)

## ----compute uppers, echo=FALSE, eval=TRUE------------------------------------
t_Rec_reconciled_mean <- c(t_rec_results$upper_rec_mean,t_rec_results$bottom_rec_mean)
t_Rec_corr_with_upper = A %*% t_rec_results$bottom_rec_scale_matrix
t_Rec_scale_par <- rbind(cbind(t_rec_results$upper_rec_scale_matrix, t_Rec_corr_with_upper),
                         cbind(t(t_Rec_corr_with_upper), t_rec_results$bottom_rec_scale_matrix))

# Index of the upper variable to plot
i_upper <- 1  # change this if a different index is needed

# Extract distribution's parameters for each method

# MinT
# Build the full reconciled covariance matrix
MinT_corr_with_upper = A %*% gauss_results$bottom_rec_cov
MinT_cov_mat = rbind(cbind(gauss_results$upper_rec_cov, MinT_corr_with_upper),
                        cbind(t(MinT_corr_with_upper),gauss_results$bottom_rec_cov))
mu_MinT <- as.numeric(MinT_reconciled_mean[i_upper])
sd_MinT <- sqrt(MinT_cov_mat[i_upper, i_upper])

# t-Rec
mu_tRec <- as.numeric(t_Rec_reconciled_mean[i_upper])
scale_tRec <- sqrt(t_Rec_scale_par[i_upper, i_upper])
df_tRec <- t_rec_results$bottom_rec_df

# Base
mu_Base <- as.numeric(Base_mean[i_upper])
sd_Base <- sqrt(Base_cov_mat[i_upper, i_upper])

# Create a grid of x values
x_vals <- seq(min(mu_MinT, mu_tRec, mu_Base) - 4*max(sd_MinT, scale_tRec, sd_Base),
              max(mu_MinT, mu_tRec, mu_Base) + 4*max(sd_MinT, scale_tRec, sd_Base),
              length.out = 1000)

# Compute densities
dens_MinT <- dnorm(x_vals, mean = mu_MinT, sd = sd_MinT)
dens_tRec <- dt((x_vals - mu_tRec) / scale_tRec, df = df_tRec) /scale_tRec
dens_Base <- dnorm(x_vals, mean = mu_Base, sd = sd_Base)


## ----comparison-plot, echo=FALSE, eval=TRUE, dpi=300, out.width = "100%", fig.align='center', fig.cap="**Figure 2**: Predictive densities of the upper time series obtained with MinT (purple), t-Rec (green) and Base (blue). The black triangle indicates the actual value.", fig.dim = c(8, 4), warning=FALSE----

dens_df <- data.frame(
  x       = rep(x_vals, 3),
  Method  = rep(c("MinT", "t-Rec", "Base"), each = length(x_vals)),
  Density = c(dens_MinT, dens_tRec, dens_Base)
)

method_colors <- c(
  "MinT" = "#440154",
  "t-Rec" = "#29AF7F",
  "Base" = "#3B528B"
)

ggplot(dens_df, aes(x = x, y = Density, color = Method)) +
  geom_area(
    data = dens_df[dens_df$Method %in% c("MinT", "t-Rec"),],
    aes(fill = Method),
    position = "identity",
    alpha = 0.3,
    color = NA
  ) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = method_colors) +
  scale_fill_manual(values = method_colors) +
  geom_point(aes(x = actuals[i_upper], y = 0, shape = "Actual value"), color = "black", size = 3) +
  scale_shape_manual(values = c("Actual value" = 17)) +
  scale_x_continuous(labels = function(x) paste0(formatC(x / 1e3, format = "g"), "k")) +
  guides(fill = "none") +
  labs(
    title = "Predictive densities of upper time series",
    x = "",
    y = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    plot.title       = element_text(size = 18)
  )

## ----compute density----------------------------------------------------------
# Select which series to plot
i = 1     # CH
j = 19    # GR

# density of the inverse gamma
dinvgamma <- function(x, shape, rate) {
  dgamma(1/x, shape = shape, rate = rate) / x^2
}

# density of the standard deviation (square root transform of variance)
d_std_dev <- function(x, shape,rate){
  dinvgamma(x^2, shape = shape, rate = rate)*2*x
}

# compute the density of the variance parameters for CH
mean_ch <- t_rec_results$posterior_Psi[i, i]/(t_rec_results$posterior_nu - n - 1)
x_ch <- sqrt(seq(mean_ch*0.5, mean_ch*1.7, length.out = 1000))
shape_ch <- (t_rec_results$posterior_nu - n + 1) /2
rate_ch  <- t_rec_results$posterior_Psi[i, i] / 2
dens_ch <- d_std_dev(x_ch, shape = shape_ch, rate = rate_ch)

# compute the density of the variance parameters for GR
mean_gr <- t_rec_results$posterior_Psi[j, j]/(t_rec_results$posterior_nu - n - 1)
x_gr <- sqrt(seq(mean_gr*0.5, mean_gr*1.7, length.out = 1000))
shape_gr <- (t_rec_results$posterior_nu - n + 1) /2
rate_gr  <- t_rec_results$posterior_Psi[j, j] / 2
dens_gr <- d_std_dev(x_gr, shape = shape_gr, rate = rate_gr)

## ----density plot, echo=FALSE, eval=TRUE, dpi=300, out.width = "100%", fig.align='center', fig.cap="**Figure 3**: Density of the posterior standard deviation of the forecasts.", fig.dim = c(8, 4), warning=FALSE----

df_dens <- data.frame(
  x = c(x_ch, x_gr), 
  y = c(dens_ch, dens_gr), 
  panel = rep(c("CH", "GR"), each = length(dens_ch))
)

# Plot the density of the standard deviation parameters for CH and GR
ggplot(df_dens, aes(x = x, y = y)) +
  geom_line(linewidth = 1.5, color="darkolivegreen", alpha=0.7) +
  scale_x_continuous(labels = function(x) paste0(formatC(x / 1e3, format = "g"), "k")) +
  guides(fill = "none") +
  labs(
    title = "Posterior standard deviation of the forecasts",
    x = "",
    y = ""
  ) +
  facet_wrap(~ panel, nrow = 1, scales = "free") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    plot.title       = element_text(size = 18)
  )

## ----generate IW samples------------------------------------------------------
# generate k samples from an IW(Psi, nu) distribution
rinvwishart <- function(k, nu, Psi, seed=42) {
  p <- nrow(Psi)
  Sigma <- solve(Psi)
  
  set.seed(seed)
  all_W <- rWishart(k, df = nu, Sigma = Sigma)
  
  W <- array(NA, dim = c(p, p, k))
  for (i in 1:k) {
    W[,,i] <- solve(all_W[,,i])
  }
  return(W)
}

IW_post_samples <- rinvwishart(k=1000, nu = t_rec_results$posterior_nu,
                            Psi = t_rec_results$posterior_Psi)


## ----compute correlations, echo=FALSE, eval=TRUE------------------------------
corr_mint <- cov2cor(MinT_cov_mat)
corr_base <- cov2cor(Base_cov_mat)
corr_tRec <- cov2cor(t_rec_results$posterior_Psi/(t_rec_results$posterior_nu - n - 1))

## ----plot densities, echo=FALSE, eval=TRUE, dpi=300, out.width = "100%", fig.align='center', fig.cap="**Figure 4**: Density of the posterior correlation between CH and GR obtained with t-Rec.", fig.dim = c(8, 4), warning=FALSE----
# Compute correlation samples for all sample matrices generated above
corr_post <- array(apply(IW_post_samples, FUN = function(M) cov2cor(M), MARGIN = c(3)),
                   dim = dim(IW_post_samples))

df_corr <- data.frame(x = corr_post[i, j, ])

ggplot(df_corr, aes(x = x)) +
  geom_density(fill = "#00BA38", alpha = 0.4) +
  geom_vline(xintercept = corr_mint[i, j], linetype = "dashed", color = "black") +
  annotate("text", x = corr_mint[i, j], y = max(density(corr_post[i, j, ])$y) * 0.95, label = "MinT", 
           hjust = -0.2, vjust = 0, size = 4.5) +
  labs(x = "", y = "Density",
       title = expression(paste("Posterior correlation ", rho["CH,GR"]))) +
  theme_minimal() +
  theme(
    legend.position  = "none",
    plot.title       = element_text(size = 18)
  )

## ----echo=FALSE,eval=TRUE-----------------------------------------------------
# Print table with values for Std and correlation
knitr::kable(data.frame(
  Method = c("Base", "MinT", "t-Rec"),
  sigma2_CH   = sqrt(c(Base_cov_mat[i, i], MinT_cov_mat[i, i], t_rec_results$posterior_Psi[i, i]/(t_rec_results$posterior_nu - n - 1))),
  sigma2_GR   = sqrt(c(Base_cov_mat[j, j], MinT_cov_mat[j, j], t_rec_results$posterior_Psi[j, j]/(t_rec_results$posterior_nu - n - 1))),
  rho_CH_GR   = c(corr_base[i, j], corr_mint[i, j], corr_tRec[i, j])
), digits = c(0, 0, 0, 2),
   format = "html",
   col.names = c("Method",
                 "$\\hat{\\sigma}_{\\text{CH}}$",
                 "$\\hat{\\sigma}_{\\text{GR}}$",
                 "$\\hat{\\rho}_{\\text{CH,GR}}$"),
   escape = FALSE,
   table.attr = "style='width:auto;'",
   format.args = list(big.mark = ",", scientific = FALSE),
   caption = "**Table 1**: Standard deviation and correlation estimates for CH and GR.")

