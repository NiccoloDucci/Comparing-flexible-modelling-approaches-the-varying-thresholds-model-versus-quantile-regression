### PREDICTION INTERVALS - model 3 normal errors i.n.i.d and model 4 normal errors i.n.i.d ###

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

#### quantile regression simulation####
coverage.rq <- function(
    x, # covariate 
    beta, # true parameter values
    err.distr, #type of error distribution
    tau, #quantiles
    xnew, # new x-value for which to predict
    nsim # number of replicates to simulate
){
  
  pi.hits <- 0
  n <- length(x)
  pred_bounds <- matrix(NA, nsim, 2)
  interval_score <- matrix(NA, nsim, 1)
  
  for (i in 1:nsim) {
    
    #error type
    err <- rnorm(n, mean = 0, sd = 1)
    
    # simulate the observed data
    if(err.distr=="het1") y <- beta[1] + beta[2]*x + exp(x-5)*err 
    if(err.distr=="het2") y <- beta[1] + beta[2]*x + exp(5-x)*err 
    
    # fit the model
    mod <- rq(formula = y~x,tau = tau )
    
    # simulate a new observation to predict
    if(err.distr=="het1"){
      ynew <- beta[1] + beta[2]*xnew + exp(xnew-5)*rnorm(1, mean = 0, sd = 1)
    }
    if(err.distr=="het2"){
      ynew <- beta[1] + beta[2]*xnew + exp(5-xnew)*rnorm(1, mean = 0, sd = 1)
    } 
    
    # compute confidence and prediction intervals
    pi <- c(mod$coefficients[1,1]+mod$coefficients[2,1]*xnew,
            mod$coefficients[1,2]+mod$coefficients[2,2]*xnew)
    if ( (pi[1] < ynew) & (ynew < pi[2])){
      pi.hits <- pi.hits + 1
    }
    pred_bounds[i,] <- pi
    interval_score[i, ] <- interval_score(
      lower_bound = pi[1],
      upper_bound = pi[2],
      observation = ynew,
      alpha = 0.05
    )
  }
  return( list( pi.hits/nsim, pred_bounds, interval_score  ) )
}

##simulation coverage 80% quantile regression
set.seed(21)
nsim=1000
nobs=1000
quantreg.coverage.model3 <-  matrix(NA, nsim, 3)
quantreg.coverage.model4 <- matrix(NA, nsim, 3)
quantreg.intScore.model3 <-  matrix(NA, nsim, 3)
quantreg.intScore.model4 <- matrix(NA, nsim, 3)

model3_bounds <- list(matrix(NA, nsim, 2), 
                      matrix(NA, nsim, 2),
                      matrix(NA, nsim, 2))

model4_bounds <-  list(matrix(NA, nsim, 2), 
                       matrix(NA, nsim, 2),
                       matrix(NA, nsim, 2))

for(k in 1:nsim){
  
  #independent variable
  mu.x <- 5
  sd.x <- 1
  x <- rnorm(nobs, mean = mu.x, sd =sd.x)
  xnew.values <- qnorm(c(.1,.5,.9), mean = mu.x, sd =sd.x)
  
  #base model
  beta0=1
  beta1=2
  
  for(i in 1:length(xnew.values)){
    
    #normal error
    model3_pred_int <- coverage.rq(x=x, 
                                   beta = c(beta0,beta1), err.distr = "het1",
                                   tau=c(0.025,0.975), xnew = xnew.values[i],
                                   nsim = 1)
    model3_bounds[[i]][k,]<- model3_pred_int[[2]]
    quantreg.coverage.model3[k,i] <- model3_pred_int[[1]]
    quantreg.intScore.model3[k,i] <- model3_pred_int[[3]]
    
    #chisq error
    model4_pred_int <- coverage.rq(x=x,
                                   beta = c(beta0,beta1), err.distr = "het2",
                                   tau=c(0.025,0.975), xnew = xnew.values[i],
                                   nsim = 1)
    model4_bounds[[i]][k,]<- model4_pred_int[[2]]
    quantreg.coverage.model4[k,i] <- model4_pred_int[[1]]
    quantreg.intScore.model4[k,i] <- model4_pred_int[[3]]
  }
  
}

# prediction interval coevrage
colMeans(quantreg.coverage.model3)
#[1] 1.000 1.000 0.923

colMeans(quantreg.coverage.model4)
# 0.908 1.000 1.000

# avg interval width
mean(model3_bounds[[1]][,2]-model3_bounds[[1]][,1])
mean(model3_bounds[[2]][,2]-model3_bounds[[2]][,1])
mean(model3_bounds[[3]][,2]-model3_bounds[[3]][,1])
# > mean(model3_bounds[[1]][,2]-model3_bounds[[1]][,1])
# [1] 3.568797
# > mean(model3_bounds[[2]][,2]-model3_bounds[[2]][,1])
# [1] 7.958519
# > mean(model3_bounds[[3]][,2]-model3_bounds[[3]][,1])
# [1] 12.30267


mean(model4_bounds[[1]][,2]-model4_bounds[[1]][,1])
mean(model4_bounds[[2]][,2]-model4_bounds[[2]][,1])
mean(model4_bounds[[3]][,2]-model4_bounds[[3]][,1])
# > mean(model4_bounds[[1]][,2]-model4_bounds[[1]][,1])
# [1] 12.2711
# > mean(model4_bounds[[2]][,2]-model4_bounds[[2]][,1])
# [1] 7.928364
# > mean(model4_bounds[[3]][,2]-model4_bounds[[3]][,1])
# [1] 3.563485

# avg interval score
colMeans(quantreg.intScore.model3)
3.568797  7.958519 16.310602

colMeans(quantreg.intScore.model4)
17.756525  7.928364  3.563485



#### varying thresholds model simulation model 3 and 4 ####

coverage.splitfit <- function(
    x, # covariate 
    beta, # true parameter values
    err.distr, #type of error distribution
    tau, #quantiles
    xnew, # new x-value for which to predict
    nsim # number of replicates to simulate
) {
  pi.hits <- 0
  n <- length(x)
  pred_bounds <- matrix(NA, nsim, 2)
  interval_score <- matrix(NA, nsim, 1)
  
  for (i in 1:nsim) {
    
    #error type
    err <- rnorm(n, mean = 0, sd = 1)
    
    # simulate the observed data
    if(err.distr=="het1") y <- beta[1] + beta[2]*x + exp(x-5)*err 
    if(err.distr=="het2") y <- beta[1] + beta[2]*x + exp(5-x)*err
    
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
    if(err.distr=="het1"){
      ynew <- beta[1] + beta[2]*xnew + exp(xnew-5)*rnorm(1, mean = 0, sd = 1)
    }
    if(err.distr=="het2"){
      ynew <- beta[1] + beta[2]*xnew + exp(5-xnew)*rnorm(1, mean = 0, sd = 1)
    } 
    
    # compute confidence and prediction intervals
    lower <- ParametricSplitsQuantiles(distribution_function = Splitfit$distribution_function,
                                       quantile = tau[1],
                                       y_values = Splitfit$interpolation_points)
    upper <- ParametricSplitsQuantiles(distribution_function = Splitfit$distribution_function,
                                       quantile = tau[2],
                                       y_values = Splitfit$interpolation_points)
    pi <- c(lower[[1]][1], upper[[1]][1])
    if ( (pi[1] < ynew) & (ynew < pi[2])){
      pi.hits <- pi.hits + 1
    }
    pred_bounds[i,] <- pi
    interval_score[i, ] <- interval_score(
      lower_bound = pi[1],
      upper_bound = pi[2],
      observation = ynew,
      alpha = 0.05
    )
  }
  return( list( pi.hits/nsim, pred_bounds, interval_score ) )
}

##simulation coverage 80% varying-threshold model 3 and model 4
set.seed(21)
nsim=1000
nobs=1000
splitfit.coverage.model3 <-  matrix(NA, nsim, 3)
splitfit.coverage.model4 <- matrix(NA, nsim, 3)
splitfit.intScore.model3 <-  matrix(NA, nsim, 3)
splitfit.intScore.model4 <- matrix(NA, nsim, 3)
splitfit.model3_bounds <- list(matrix(NA, nsim, 2), 
                               matrix(NA, nsim, 2),
                               matrix(NA, nsim, 2))

splitfit.model4_bounds <-  list(matrix(NA, nsim, 2), 
                                matrix(NA, nsim, 2),
                                matrix(NA, nsim, 2))

for(k in 1:nsim){
  
  #independent variable
  mu.x <- 5
  sd.x <- 1
  x <- rnorm(nobs, mean = mu.x, sd =sd.x)
  xnew.values <- qnorm(c(.1,.5,.9), mean = mu.x, sd =sd.x)
  
  #base model
  beta0=1
  beta1=2
  
  for(i in 1:length(xnew.values)){
    
    #het1 error
    splitfit.model3_pred_int <- coverage.splitfit(x=x, 
                                                  beta = c(beta0,beta1), 
                                                  err.distr = "het1",
                                                  tau=c(0.025,0.975), 
                                                  xnew = xnew.values[i],
                                                  nsim = 1)
    
    splitfit.model3_bounds[[i]][k,] <-  splitfit.model3_pred_int[[2]]
    splitfit.coverage.model3[k,i] <- splitfit.model3_pred_int[[1]]
    splitfit.intScore.model3[k,i] <- splitfit.model3_pred_int[[3]]
    
    #het2 error
    splitfit.model4_pred_int <- coverage.splitfit(x=x,
                                                  beta = c(beta0,beta1), 
                                                  err.distr = "het2",
                                                  tau=c(0.025,0.975), 
                                                  xnew = xnew.values[i],
                                                  nsim = 1)
    splitfit.coverage.model4[k,i] <- splitfit.model4_pred_int[[1]]
    splitfit.model4_bounds[[i]][k,] <- splitfit.model4_pred_int[[2]]
    splitfit.intScore.model4[k,i] <- splitfit.model4_pred_int[[3]]
  }
  print(k)
}


# # coverage
# colMeans(splitfit.coverage.model3, na.rm = F)
# [1] 1.000 0.998 0.899
# 
# colMeans(splitfit.coverage.model4, na.rm = F )
# [1] 0.902 0.999 1.000
# 
# # interval width
# mean( splitfit.model4_bounds[[1]][,2] - splitfit.model4_bounds[[1]][,1])
# [1] 18.47979
# mean( splitfit.model4_bounds[[2]][,2] - splitfit.model4_bounds[[2]][,1])
# [1] 7.400119
# mean( splitfit.model4_bounds[[3]][,2] - splitfit.model4_bounds[[3]][,1])
# [1] 6.822311
# 
# 
# mean( splitfit.model3_bounds[[1]][,2] - splitfit.model3_bounds[[1]][,1])
# [1] 6.669608
# mean( splitfit.model3_bounds[[2]][,2] - splitfit.model3_bounds[[2]][,1])
# [1] 7.391654
# mean( splitfit.model3_bounds[[3]][,2] - splitfit.model3_bounds[[3]][,1])
# [1] 18.4868
# 
# # interval score
# colMeans(splitfit.intScore.model3)
# [1]  6.669608  7.411906 24.402656
# 
# 
# colMeans(splitfit.intScore.model4)
# [1] 24.914052  7.428959  6.822311



#####

