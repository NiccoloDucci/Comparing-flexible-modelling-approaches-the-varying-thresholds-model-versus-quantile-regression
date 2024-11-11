### PREDICTION INTERVALS - model 6 cubic effect ###

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

#### quantile regression simulation ####

coverage.rq.cubic <- function(
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
    # simulate the observed data
    y <- beta[1] + beta[2] * (x) + beta[3] * (x^2) + beta[4] * (x^3) + err
    # fit the model
    mod <- rq(formula = y ~ x, tau = tau)
    # simulate a new observation to predict
    if (err.distr == "normal") {
      ynew <- beta[1] + beta[2] * (xnew) + beta[3] * (xnew^2) + beta[4] * (xnew^3) + rnorm(1, 0, 1)
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
      alpha = 0.05
    )
  }
  return(list(pi.hits / nsim, pred_bounds, interval_score))
}

## simulation coverage 80% quantile regression i.i.d case
set.seed(21)
nsim <- 500
nobs <- 1000
quantreg.coverage.model6 <- matrix(NA, nsim, 3)
quantreg.intScore.model6 <- matrix(NA, nsim, 3)
model6_bounds <- list(
  matrix(NA, nsim, 2),
  matrix(NA, nsim, 2),
  matrix(NA, nsim, 2)
)

for (k in 1:nsim) {
  # independent variable
  max.x <- 8
  min.x <- -3
  x <- runif(nobs, min = min.x, max = max.x)
  xnew.values <- qunif(c(.1, .5, .9), min = min.x, max = max.x)

  # base model
  beta <- c(-0.2, 1, -0.8, 0.12)

  for (i in 1:length(xnew.values)) {
    # normal error
    model6_pred_int <- coverage.rq.cubic(
      x = x,
      beta = beta, err.distr = "normal",
      tau = c(0.05, 0.975), xnew = xnew.values[i],
      nsim = 1
    )
    model6_bounds[[i]][k, ] <- model6_pred_int[[2]]
    quantreg.coverage.model6[k, i] <- model6_pred_int[[1]]
    quantreg.intScore.model6[k, i] <- model6_pred_int[[3]]
  }
}

# coverage
colMeans(quantreg.coverage.model6)

# avg interval width
mean(model6_bounds[[1]][, 2] - model6_bounds[[1]][, 1])
mean(model6_bounds[[2]][, 2] - model6_bounds[[2]][, 1])
mean(model6_bounds[[3]][, 2] - model6_bounds[[3]][, 1])

# avg interval score
colMeans(quantreg.intScore.model6)



#### varying thresholds model simulation model 6 ####

coverage.splitfit.cubic <- function(
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

    # simulate the observed data
    y <- beta[1] + beta[2] * (x) + beta[3] * (x^2) + beta[4] * (x^3) + err

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

    # VTM
    Splitfit <- VTM(
      data_train = ds,
      data_predict = datp,
      formula = formbin,
      splits = splits,
      num_interpolation_points = length(splits) * 2,
      family_link = binomial(link = "probit"),
      lambda = 0,
      alpha = 0
    )

    # simulate a new observation to predict
    if (err.distr == "normal") {
      ynew <- beta[1] + beta[2] * (xnew) + beta[3] * (xnew^2) + beta[4] * (xnew^3) + rnorm(1, mean = 0, sd = 1)
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
      alpha = 0.05
    )
  }
  return(list(pi.hits / nsim, pred_bounds, interval_score))
}

## simulation coverage 80% varying-threshold model i.i.d case
set.seed(21)
nsim <- 500
nobs <- 1000
splitfit.coverage.model6 <- matrix(NA, nsim, 3)
splitfit.intScore.model6 <- matrix(NA, nsim, 3)
splitfit.model6_bounds <- list(
  matrix(NA, nsim, 2),
  matrix(NA, nsim, 2),
  matrix(NA, nsim, 2)
)


for (k in 1:nsim) {
  # independent variable
  max.x <- 8
  min.x <- -3
  x <- runif(nobs, min = min.x, max = max.x)
  xnew.values <- qunif(c(.1, .5, .9), min = min.x, max = max.x)

  # base model
  beta <- c(-0.2, 1, -0.8, 0.12)

  for (i in 1:length(xnew.values)) {
    # normal error
    splitfit.model6_pred_int <- coverage.splitfit.cubic(
      x = x,
      beta = beta,
      err.distr = "normal",
      tau = c(0.05, 0.975),
      xnew = xnew.values[i],
      nsim = 1
    )

    splitfit.model6_bounds[[i]][k, ] <- splitfit.model6_pred_int[[2]]
    splitfit.coverage.model6[k, i] <- splitfit.model6_pred_int[[1]]
    splitfit.intScore.model6[k, i] <- splitfit.model6_pred_int[[3]]
  }
  print(k)
}

# coverage
colMeans(splitfit.coverage.model6, na.rm = F)


# avg interval width
mean(splitfit.model6_bounds[[1]][, 2] - splitfit.model6_bounds[[1]][, 1])
mean(splitfit.model6_bounds[[2]][, 2] - splitfit.model6_bounds[[2]][, 1])
mean(splitfit.model6_bounds[[3]][, 2] - splitfit.model6_bounds[[3]][, 1])


# avg interval score
colMeans(splitfit.intScore.model6, na.rm = F)
