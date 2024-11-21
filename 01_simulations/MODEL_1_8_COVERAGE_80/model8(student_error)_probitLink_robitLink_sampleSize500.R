### PREDICTION INTERVALS - model 8 t student (df=3) i.i.d. errors with quadratic effect ###

##### Libraries#####
library(tidyverse)
library(ggplot2)
library(quantreg)
library(VGAM)
library(glmx)
library(scoringRules)
#####

##### functions #####

source("C:\\Users\\nicco\\Desktop\\UNI\\00Esami in corso\\Tesi\\ADAC_revision\\00_functions\\varying_thresholds_model.R")
source("C:\\Users\\nicco\\Desktop\\UNI\\00Esami in corso\\Tesi\\ADAC_revision\\00_functions\\parametric_splits_quantiles.R")
source("C:\\Users\\nicco\\Desktop\\UNI\\00Esami in corso\\Tesi\\ADAC_revision\\00_functions\\interval_score.R")

#####

## PREDICTION INTERVAL SIMULATION ##

#### quantile regression simulation ####

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
    if (err.distr == "student") {
      err <- rt(n, df = 3)
    }
    # simulate the observed data
    y <- beta[1] + beta[2] * x + beta[3] * (x^2) + err
    # fit the model
    mod <- rq(formula = y ~ x, tau = tau)
    # simulate a new observation to predict
    if (err.distr == "student") {
      ynew <- beta[1] + beta[2] * xnew + beta[3] * (xnew^2) + rt(1, df = 3)
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
nobs <- 500
quantreg.coverage.model8 <- matrix(NA, nsim, 3)
quantreg.intScore.model8 <- matrix(NA, nsim, 3)
model8_bounds <- list(
  matrix(NA, nsim, 2),
  matrix(NA, nsim, 2),
  matrix(NA, nsim, 2)
)

for (k in 1:nsim) {
  # independent variable
  x <- runif(nobs, 0, 6)
  xnew.values <- qunif(c(.1, .5, .9), 0, 6)
  
  # base model
  beta0 <- 1
  beta1 <- 0.5
  beta2 <- 1
  
  for (i in 1:length(xnew.values)) {
    # student error(df=3)
    model8_pred_int <- coverage.rq(
      x = x,
      beta = c(beta0, beta1, beta2), err.distr = "student",
      tau = c(0.1, 0.9), xnew = xnew.values[i],
      nsim = 1
    )
    model8_bounds[[i]][k, ] <- model8_pred_int[[2]]
    quantreg.coverage.model8[k, i] <- model8_pred_int[[1]]
    quantreg.intScore.model8[k, i] <- model8_pred_int[[3]]
  }
}

# coverage
colMeans(quantreg.coverage.model8)


# avg interval width
mean(model8_bounds[[1]][, 2] - model8_bounds[[1]][, 1])
mean(model8_bounds[[2]][, 2] - model8_bounds[[2]][, 1])
mean(model8_bounds[[3]][, 2] - model8_bounds[[3]][, 1])


# avg interval score
colMeans(quantreg.intScore.model8)


#### PROBIT LINK: varying thresholds model simulation model 8  ####

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
    if (err.distr == "student") {
      err <- rt(n, df = 3)
    }
    # simulate the observed data
    y <- beta[1] + beta[2] * x + beta[3] * (x^2) + err
    
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
    if (err.distr == "student") {
      ynew <- beta[1] + beta[2] * xnew + beta[3] * (xnew^2) + rt(1, df = 3)
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
nobs <- 500
splitfit.coverage.model8 <- matrix(NA, nsim, 3)
splitfit.intScore.model8 <- matrix(NA, nsim, 3)
splitfit.model8_bounds <- list(
  matrix(NA, nsim, 2),
  matrix(NA, nsim, 2),
  matrix(NA, nsim, 2)
)

for (k in 1:nsim) {
  # independent variable
  x <- runif(nobs, 0, 6)
  xnew.values <- qunif(c(.1, .5, .9), 0, 6)
  
  # base model
  beta0 <- 1
  beta1 <- 0.5
  beta2 <- 1
  
  for (i in 1:length(xnew.values)) {
    # student error(df=3)
    splitfit.model8_pred_int <- coverage.splitfit(
      x = x,
      beta = c(beta0, beta1, beta2),
      err.distr = "student",
      tau = c(0.1, 0.9),
      xnew = xnew.values[i],
      nsim = 1
    )
    splitfit.coverage.model8[k, i] <- splitfit.model8_pred_int[[1]]
    splitfit.model8_bounds[[i]][k, ] <- splitfit.model8_pred_int[[2]]
    splitfit.intScore.model8[k, i] <- splitfit.model8_pred_int[[3]]
  }
  print(k)
}


> # coverage
  > colMeans(splitfit.coverage.model8)
[1] 0.786 0.874 0.841
> 
  > 
  > # avg interval width
  > mean(splitfit.model8_bounds[[1]][, 2] - splitfit.model8_bounds[[1]][, 1])
[1] 3.359758
> mean(splitfit.model8_bounds[[2]][, 2] - splitfit.model8_bounds[[2]][, 1])
[1] 4.936328
> mean(splitfit.model8_bounds[[3]][, 2] - splitfit.model8_bounds[[3]][, 1])
[1] 4.425607
> 
  > # avg interval score
  > colMeans(splitfit.intScore.model8)
[1] 5.820284 6.961117 6.288258



#### ROBIT LINK: varying thresholds model simulation model 8  ####

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
    if (err.distr == "student") {
      err <- rt(n, df = 3)
    }
    # simulate the observed data
    y <- beta[1] + beta[2] * x + beta[3] * (x^2) + err
    
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
    # robit link
    link_new <- gosset(nu = 3)
    
    # VTM
    Splitfit <- VTM(
      data_train = ds,
      data_predict = datp,
      formula = formbin,
      splits = splits,
      num_interpolation_points = length(splits) * 2,
      family_link = binomial(link = link_new),
      lambda = 0,
      alpha = 0
    )
    
    
    # simulate a new observation to predict
    if (err.distr == "student") {
      ynew <- beta[1] + beta[2] * xnew + beta[3] * (xnew^2) + rt(1, df = 3)
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
nobs <- 500
splitfit.coverage.model8_probit <- matrix(NA, nsim, 3)
splitfit.intScore.model8_probit <- matrix(NA, nsim, 3)
splitfit.model8_bounds_probit <- list(
  matrix(NA, nsim, 2),
  matrix(NA, nsim, 2),
  matrix(NA, nsim, 2)
)

for (k in 1:nsim) {
  # independent variable
  x <- runif(nobs, 0, 6)
  xnew.values <- qunif(c(.1, .5, .9), 0, 6)
  
  # base model
  beta0 <- 1
  beta1 <- 0.5
  beta2 <- 1
  
  for (i in 1:length(xnew.values)) {
    # student error(df=3)
    splitfit.model8_pred_int_probit <- coverage.splitfit(
      x = x,
      beta = c(beta0, beta1, beta2),
      err.distr = "student",
      tau = c(0.1, 0.9),
      xnew = xnew.values[i],
      nsim = 1
    )
    splitfit.coverage.model8_probit[k, i] <- splitfit.model8_pred_int_probit[[1]]
    splitfit.model8_bounds_probit[[i]][k, ] <- splitfit.model8_pred_int_probit[[2]]
    splitfit.intScore.model8_probit[k, i] <- splitfit.model8_pred_int_probit[[3]]
  }
  print(k)
}

> # coverage
  > colMeans(splitfit.coverage.model8_probit, na.rm = F)
[1] 0.745 0.773 0.757
> 
  > # avg interval width
  > mean(splitfit.model8_bounds_probit[[1]][, 2] - splitfit.model8_bounds_probit[[1]][, 1])
[1] 2.938146
> mean(splitfit.model8_bounds_probit[[2]][, 2] - splitfit.model8_bounds_probit[[2]][, 1])
[1] 3.367132
> mean(splitfit.model8_bounds_probit[[3]][, 2] - splitfit.model8_bounds_probit[[3]][, 1])
[1] 3.306127
> 
  > # avg interval score
  > colMeans(splitfit.intScore.model8_probit, na.rm = F)
[1] 5.894800 6.461548 5.870803
