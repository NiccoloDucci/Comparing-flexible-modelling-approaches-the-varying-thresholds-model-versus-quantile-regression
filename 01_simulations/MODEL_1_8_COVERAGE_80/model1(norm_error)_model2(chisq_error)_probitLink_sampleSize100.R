### PREDICTION INTERVALS - model 1 normal errors i.i.d and model 2 chi-squared errors i.i.d ###

##### Libraries #####
library(tidyverse)
library(ggplot2)
library(quantreg)
#####

##### functions #####

source("C:\\Users\\nicco\\Desktop\\UNI\\00Esami in corso\\Tesi\\ADAC_revision\\00_functions\\varying_thresholds_model.R")
source("C:\\Users\\nicco\\Desktop\\UNI\\00Esami in corso\\Tesi\\ADAC_revision\\00_functions\\parametric_splits_quantiles.R")
source("C:\\Users\\nicco\\Desktop\\UNI\\00Esami in corso\\Tesi\\ADAC_revision\\00_functions\\interval_score.R")

#####


## PREDICTION INTERVAL SIMULATION ##

#### quantile regression  simulation####
coverage.rq <- function(
    x, # covariate
    beta, # true parameter values
    err.distr, # type of error distribution
    tau, # quantiles
    xnew, # new x-value for which to predict
    nsim # number of replicates to simulate
) {
  pi.hits <- 0
  n <- length(x)
  pred_bounds <- matrix(NA, nsim, 2)
  interval_score <- matrix(NA, nsim, 1)
  for (i in 1:nsim) {
    # error type
    if (err.distr == "normal") {
      err <- rnorm(n, mean = 0, sd = 1)
    }
    if (err.distr == "chisq") {
      err <- rchisq(n, df = 3)
    }
    # simulate the observed data
    y <- beta[1] + beta[2] * x + err
    # fit the model
    mod <- rq(formula = y ~ x, tau = tau)
    # simulate a new observation to predict
    if (err.distr == "normal") {
      ynew <- beta[1] + beta[2] * xnew + rnorm(1, mean = 0, sd = 1)
    }
    if (err.distr == "chisq") {
      ynew <- beta[1] + beta[2] * xnew + rchisq(1, df = 3)
    }
    # compute confidence and prediction intervals
    pi <- c(
      mod$coefficients[1, 1] + mod$coefficients[2, 1] * xnew,
      mod$coefficients[1, 2] + mod$coefficients[2, 2] * xnew
    )
    if ((pi[1] < ynew) & (ynew < pi[2])) {
      pi.hits <- pi.hits + 1
    }
    pred_bounds[i, ] <- pi
    interval_score[i, ] <- interval_score(
      lower_bound = pi[1],
      upper_bound = pi[2],
      observation = ynew,
      alpha = 0.2
    )
  }
  return(list(pi.hits / nsim, pred_bounds, interval_score))
}

## simulation coverage 80% quantile regression i.i.d case
set.seed(21)
nsim <- 1000
nobs <- 100
quantreg.coverage.norm <- matrix(NA, nsim, 3)
quantreg.intScore.norm <- matrix(NA, nsim, 3)
quantreg.coverage.chisq <- matrix(NA, nsim, 3)
quantreg.intScore.chisq <- matrix(NA, nsim, 3)

normal_bounds <- list(
  matrix(NA, nsim, 2),
  matrix(NA, nsim, 2),
  matrix(NA, nsim, 2)
)

chisq_bounds <- list(
  matrix(NA, nsim, 2),
  matrix(NA, nsim, 2),
  matrix(NA, nsim, 2)
)

for (k in 1:nsim) {
  # independent variable
  mu.x <- 5
  sd.x <- 1
  x <- rnorm(nobs, mean = mu.x, sd = sd.x)
  xnew.values <- qnorm(c(.1, .5, .9), mean = mu.x, sd = sd.x)
  
  # base model
  beta0 <- 1
  beta1 <- 2
  
  for (i in 1:length(xnew.values)) {
    # normal error
    normal_pred_int <- coverage.rq(
      x = x,
      beta = c(beta0, beta1), err.distr = "normal",
      tau = c(0.1, 0.9), xnew = xnew.values[i],
      nsim = 1
    )
    normal_bounds[[i]][k, ] <- normal_pred_int[[2]]
    quantreg.coverage.norm[k, i] <- normal_pred_int[[1]]
    quantreg.intScore.norm[k, i] <- normal_pred_int[[3]]
    
    # chisq error
    chisq_pred_int <- coverage.rq(
      x = x,
      beta = c(beta0, beta1), err.distr = "chisq",
      tau = c(0.1, 0.9), xnew = xnew.values[i],
      nsim = 1
    )
    chisq_bounds[[i]][k, ] <- chisq_pred_int[[2]]
    quantreg.coverage.chisq[k, i] <- chisq_pred_int[[1]]
    quantreg.intScore.chisq[k, i] <- chisq_pred_int[[3]]
  }
}

# model 1 coverage
colMeans(quantreg.coverage.norm)
# 0.774 0.791 0.783

# model 2 coverage
colMeans(quantreg.coverage.chisq)
# 0.793 0.746 0.765



# model 1 normal interval width
mean(normal_bounds[[1]][, 2] - normal_bounds[[1]][, 1])
mean(normal_bounds[[2]][, 2] - normal_bounds[[2]][, 1])
mean(normal_bounds[[3]][, 2] - normal_bounds[[3]][, 1])
#> mean(normal_bounds[[1]][, 2] - normal_bounds[[1]][, 1])
#[1] 2.547
#> mean(normal_bounds[[2]][, 2] - normal_bounds[[2]][, 1])
#[1] 2.557
#> mean(normal_bounds[[3]][, 2] - normal_bounds[[3]][, 1])
#[1] 2.543

# model 2 chi-square interval width
mean(chisq_bounds[[1]][, 2] - chisq_bounds[[1]][, 1])
mean(chisq_bounds[[2]][, 2] - chisq_bounds[[2]][, 1])
mean(chisq_bounds[[3]][, 2] - chisq_bounds[[3]][, 1])
#> mean(chisq_bounds[[1]][, 2] - chisq_bounds[[1]][, 1])
#[1] 5.640
#> mean(chisq_bounds[[2]][, 2] - chisq_bounds[[2]][, 1])
#[1] 5.674
#> mean(chisq_bounds[[3]][, 2] - chisq_bounds[[3]][, 1])
#[1] 5.634

# model 1 normal interval score
colMeans(quantreg.intScore.norm)
# 3.661699 3.547086 3.432073

# model 2 chi-square interval score
colMeans(quantreg.intScore.chisq)
# 8.450921 8.631527 8.702347

#####

#### varying threshold model simulation####

coverage.splitfit <- function(
    x, # covariate
    beta, # true parameter values
    err.distr, # type of error distribution
    tau, # quantiles
    xnew, # new x-value for which to predict
    nsim # number of replicates to simulate
) {
  pi.hits <- 0
  n <- length(x)
  pred_bounds <- matrix(NA, nsim, 2)
  interval_score <- matrix(NA, nsim, 1)
  
  for (i in 1:nsim) {
    # error type
    if (err.distr == "normal") {
      err <- rnorm(n, mean = 0, sd = 1)
    }
    if (err.distr == "chisq") {
      err <- rchisq(n, df = 3)
    }
    # simulate the observed data
    y <- beta[1] + beta[2] * x + err
    
    # fit the model
    ### splits generation
    # k=52 (50 + 2), there are k+1 theta thresholds
    splitnumber <- 50
    # theta1 - thetak-1
    minsplit <- quantile(y, probs = 0.02)
    maxsplit <- quantile(y, probs = 0.98)
    increase <- (maxsplit - minsplit) / splitnumber
    splits <- seq(minsplit, maxsplit, by = increase)
    # formula
    formbin <- binresp ~ x
    # dataset modification
    ds <- data.frame(y = y, x = x)
    names(ds)[1] <- "resp"
    datp <- data.frame(x = xnew) ### sample to predict
    
    #VTM
    Splitfit <- VTM(
      data_train = ds,
      data_predict = datp,
      formula = formbin,
      splits = splits,
      num_interpolation_points = length(splits) * 2,
      family_link = binomial( link = "probit" ),
      lambda = 0,
      alpha = 0
    )
    
    # simulate a new observation to predict
    if (err.distr == "normal") {
      ynew <- beta[1] + beta[2] * xnew + rnorm(1, mean = 0, sd = 1)
    }
    if (err.distr == "chisq") {
      ynew <- beta[1] + beta[2] * xnew + rchisq(1, df = 3)
    }
    # compute confidence and prediction intervals
    lower <- ParametricSplitsQuantiles(
      distribution_function = Splitfit$distribution_function,
      quantile = tau[1], y_values = Splitfit$interpolation_points
    )
    upper <- ParametricSplitsQuantiles(
      distribution_function = Splitfit$distribution_function,
      quantile = tau[2], y_values = Splitfit$interpolation_points
    )
    pi <- c(lower[[1]][1], upper[[1]][1])
    if ((pi[1] < ynew) & (ynew < pi[2])) {
      pi.hits <- pi.hits + 1
    }
    pred_bounds[i, ] <- pi
    interval_score[i, ] <- interval_score(
      lower_bound = pi[1],
      upper_bound = pi[2],
      observation = ynew,
      alpha = 0.2
    )
  }
  return(list(pi.hits / nsim, pred_bounds, interval_score))
}

## simulation coverage 80% varying-threshold model i.i.d case
set.seed(21)
nsim <- 1000
nobs <- 100
splitfit.coverage.norm <- matrix(NA, nsim, 3)
splitfit.coverage.chisq <- matrix(NA, nsim, 3)

splitfit.intScore.norm <- matrix(NA, nsim, 3)
splitfit.intScore.chisq <- matrix(NA, nsim, 3)

splitfit.normal_bounds <- list(
  matrix(NA, nsim, 2),
  matrix(NA, nsim, 2),
  matrix(NA, nsim, 2)
)
splitfit.chisq_bounds <- list(
  matrix(NA, nsim, 2),
  matrix(NA, nsim, 2),
  matrix(NA, nsim, 2)
)

for (k in 1:nsim) {
  # independent variable
  mu.x <- 5
  sd.x <- 1
  x <- rnorm(nobs, mean = mu.x, sd = sd.x)
  xnew.values <- qnorm(c(.1, .5, .9), mean = mu.x, sd = sd.x)
  
  # base model
  beta0 <- 1
  beta1 <- 2
  
  for (i in 1:length(xnew.values)) {
    # normal error
    splitfit.normal_pred_int <- coverage.splitfit(
      x = x,
      beta = c(beta0, beta1),
      err.distr = "normal",
      tau = c(0.1, 0.9),
      xnew = xnew.values[i],
      nsim = 1
    )
    
    splitfit.normal_bounds[[i]][k, ] <- splitfit.normal_pred_int[[2]]
    splitfit.coverage.norm[k, i] <- splitfit.normal_pred_int[[1]]
    splitfit.intScore.norm[k, i] <- splitfit.normal_pred_int[[3]]
    
    # chisq error
    splitfit.chisq_pred_int <- coverage.splitfit(
      x = x,
      beta = c(beta0, beta1),
      err.distr = "chisq",
      tau = c(0.1, 0.9),
      xnew = xnew.values[i],
      nsim = 1
    )
    splitfit.coverage.chisq[k, i] <- splitfit.chisq_pred_int[[1]]
    splitfit.chisq_bounds[[i]][k, ] <- splitfit.chisq_pred_int[[2]]
    splitfit.intScore.chisq[k, i] <- splitfit.chisq_pred_int[[3]]
  }
  print(k)
}

# model 1 coverage
colMeans(splitfit.coverage.norm)
#  0.812 0.773 0.788

# model 2 coverage
colMeans(splitfit.coverage.chisq)
# 0.799 0.739 0.882

# model 1 normal interval width
mean(splitfit.normal_bounds[[1]][, 2] - splitfit.normal_bounds[[1]][, 1])
mean(splitfit.normal_bounds[[2]][, 2] - splitfit.normal_bounds[[2]][, 1])
mean(splitfit.normal_bounds[[3]][, 2] - splitfit.normal_bounds[[3]][, 1])
#> mean(splitfit.normal_bounds[[1]][, 2] - splitfit.normal_bounds[[1]][, 1])
#[1] 2.98312
#> mean(splitfit.normal_bounds[[2]][, 2] - splitfit.normal_bounds[[2]][, 1])
#[1] 2.475637
#> mean(splitfit.normal_bounds[[3]][, 2] - splitfit.normal_bounds[[3]][, 1])
#[1] 2.96434


# model 2 chi-square interval width
mean(splitfit.chisq_bounds[[1]][, 2] - splitfit.chisq_bounds[[1]][, 1])
mean(splitfit.chisq_bounds[[2]][, 2] - splitfit.chisq_bounds[[2]][, 1])
mean(splitfit.chisq_bounds[[3]][, 2] - splitfit.chisq_bounds[[3]][, 1])
#> mean(splitfit.chisq_bounds[[1]][, 2] - splitfit.chisq_bounds[[1]][, 1])
#[1] 5.300
#> mean(splitfit.chisq_bounds[[2]][, 2] - splitfit.chisq_bounds[[2]][, 1])
#[1] 5.688
#> mean(splitfit.chisq_bounds[[3]][, 2] - splitfit.chisq_bounds[[3]][, 1])
#[1] 7.451

# model 1 normal interval score
colMeans(splitfit.intScore.norm)
#  4.098365 3.663728 3.977511

# model 2 chi-square interval score
colMeans(splitfit.intScore.chisq)
# 9.026443 8.713211 9.317977


#####
