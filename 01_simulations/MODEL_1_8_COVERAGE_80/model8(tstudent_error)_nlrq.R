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



#### POLYNOMIAL QUANTILE REGRESSION degree 2 - model 8####


# Polynomial regression (2nd order)
poly2.fun <- function(predictor, a, b, c) {
  a + b * predictor + c * (predictor^2)
}

poly2.Init <- function(mCall, LHS, data, ...) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  lmFit <- lm((xy[, "y"]) ~ xy[, "x"] + I(xy[, "x"]^2))
  coefs <- coef(lmFit)
  a <- coefs[1]
  b <- coefs[2]
  c <- coefs[3]
  value <- c(a, b, c)
  names(value) <- mCall[c("a", "b", "c")]
  value
}
NLS.poly2 <- selfStart(poly2.fun, poly2.Init, parameters = c("a", "b", "c"))



## quantile regression
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
    if (err.distr == "student") {
      err <- rt(n, df = 3)
    }
    # simulate the observed data
    y <- beta[1] + beta[2] * x + beta[3] * (x^2) + err

    # fit the model
    mod1 <- nlrq(
      formula = y ~ NLS.poly2(x, a, b, c),
      tau = tau[1],
      trace = T
    )

    mod2 <- nlrq(
      formula = y ~ NLS.poly2(x, a, b, c),
      tau = tau[2],
      trace = T
    )
    # simulate a new observation to predict
    if (err.distr == "student") {
      ynew <- beta[1] + beta[2] * xnew + beta[3] * (xnew^2) + rt(1, df = 3)
    }
    # compute confidence and prediction intervals
    pi <- c(
      mod1$m$getPars() %*% c(1, xnew, xnew^2),
      mod2$m$getPars() %*% c(1, xnew, xnew^2)
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
nobs <- 1000
quantreg.coverage.poly2 <- matrix(NA, nsim, 3)
quantreg.intScore.poly2 <- matrix(NA, nsim, 3)
poly2_bounds <- list(
  matrix(NA, nsim, 2),
  matrix(NA, nsim, 2),
  matrix(NA, nsim, 2)
)

for (k in 1:nsim) {
  
  #independent variable
  x <- runif(nobs,0,6)
  xnew.values <- qunif(c(.1,.5,.9), 0,6)
  
  #base model
  beta0=1
  beta1=0.5
  beta2=1
  beta = c(beta0, beta1, beta2)
  
  for (i in 1:length(xnew.values)) {
    # student error(df=3)
    poly2_pred_int <- coverage.rq.cubic(
      x = x,
      beta = beta, err.distr = "student",
      tau = c(0.1, 0.9), xnew = xnew.values[i],
      nsim = 1
    )
    poly2_bounds[[i]][k, ] <- poly2_pred_int[[2]]
    quantreg.coverage.poly2[k, i] <- poly2_pred_int[[1]]
    quantreg.intScore.poly2[k, i] <- poly2_pred_int[[3]]
  }
}

# coverage
colMeans(quantreg.coverage.poly2)
# 0.786 0.800 0.806

# interval width
mean(poly2_bounds[[1]][, 2] - poly2_bounds[[1]][, 1])
mean(poly2_bounds[[2]][, 2] - poly2_bounds[[2]][, 1])
mean(poly2_bounds[[3]][, 2] - poly2_bounds[[3]][, 1])
# mean(poly2_bounds[[1]][, 2] - poly2_bounds[[1]][, 1])
#[1] 3.265423
# mean(poly2_bounds[[2]][, 2] - poly2_bounds[[2]][, 1])
#[1] 3.268109
# mean(poly2_bounds[[3]][, 2] - poly2_bounds[[3]][, 1])
#[1] 3.269678


# interval score
colMeans(quantreg.intScore.poly2)
# 5.801 5.524 5.457

#####
