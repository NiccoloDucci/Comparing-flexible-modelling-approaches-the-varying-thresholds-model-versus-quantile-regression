###PREDICTION INTERVALS- non-parametric regression model 8###

##### functions #####

source("C:\\Users\\nicco\\Desktop\\UNI\\00Esami in corso\\Tesi\\ADAC_revision\\00_functions\\varying_thresholds_model.R")
source("C:\\Users\\nicco\\Desktop\\UNI\\00Esami in corso\\Tesi\\ADAC_revision\\00_functions\\parametric_splits_quantiles.R")
source("C:\\Users\\nicco\\Desktop\\UNI\\00Esami in corso\\Tesi\\ADAC_revision\\00_functions\\interval_score.R")

#####

#####Libraries#####
library(tidyverse)
library(ggplot2)
library(quantreg)
library(VGAM)
#####

#####SIMULATED DATASETS - model 8#####

#seed
set.seed(21)

#independent variable
max.x <- 8
min.x <- -3
x <- runif(nobs,min = min.x, max = max.x )

#errors standard normal
err_normal <- rnorm(nobs, 0, 1)

#base model
beta0=-0.2
beta1=1
beta2=-0.8
beta3=0.12
y=beta0 + beta1*(x) + beta2*(x^2) + beta3*(x^3)

#normal model
y.norm <- y+err_normal

#cubic
datset.cubic <- data.frame(y=y.norm, x=x)

#####

####Choose best lambda on a random sample#####

MCV <- function(lambdas, formula, data, tau = 0.5, k = 5){
  F <- Munge(formula, lambdas = lambdas)
  f <- rqss(F, data, tau = tau)
  n <- f$n
  m <- length(f$qss)
  y <- f$y[1:n]
  folds = sample(rep(1:k, length = n))
  U = NULL
  for(i in 1:k){
    s = which(folds != i)
    M = rqss(F, data = data[s,], tau = tau)
    nd = data[-s,]
    G = matrix(0,nrow(nd),m)
    for(j in 1:m){ #remove extrapolates, if any
      g = data$x
      G[,j] = (min(g[s]) < g[-s]) & (g[-s] < max(g[s]))
    }
    h = as.logical(apply(G,1,prod))
    u = predict(M, newdata = nd[h,]) - (y[-s])[h]
    U = c(U,(u * (tau - (u < 0))))
  }
  mean(U)
}


set.seed(22)
lams <- mcvs <- seq(0.1, 2.5, by = 0.2)
for(i in 1:length(mcvs))
  mcvs[i] <- MCV(lams[i], y ~ qss(x, lambda = lambdas[1]),
                 datset.norm)
par(mfrow = c(1,2))
plot(lams, mcvs, cex = .5, lwd = 2, type = 'l',
     xlab = expression(lambda), ylab = expression(MCV( lambda )))
lambdastar <- lams[which.min(mcvs)]

plot(datset.norm$x, datset.norm$y, cex = .5, col = "grey")
f <- rqss(y ~ qss(x, lambda = lambdastar),
          data = datset.norm, tau=0.5)
plot(f, add = TRUE, lwd = 2)
text(10, 1,bquote(lambda == ~ .(lambdastar)))

#lambda
lambda=lambdastar

#rq splines
mod1 <- rqss(y ~ qss(x, constraint= "N",lambda = lambda), 
             data = datset.norm,
             tau=0.1)

mod2 <- rqss(y ~ qss(x, constraint= "N",lambda = lambda), 
             data = datset.norm,
             tau=0.9)
ggplot()+
  geom_point(aes(x=datset.cubic$x, y=datset.cubic$y))+
  geom_line(aes(x=datset.cubic$x, 
                y= predict.rqss(mod1,datset.cubic)),
            col="red",
            lwd=1.6)+
  geom_line(aes(x=datset.cubic$x, 
                y=predict.rqss(mod2,datset.cubic)),
            col="red",
            lwd=1.6)+
  theme_classic()
#####


#####GOSSET PREDICTION INTERVAL SIMULATION - cubic #####

##quantile regression
coverage.rq.cubic <- function(
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
    if(err.distr=="normal"){
      err <- rnorm(n, mean = 0, sd = 1)
    }
    # simulate the observed data
    y=beta[1] + beta[2]*(x) + beta[3]*(x^2) + beta[4]*(x^3) + err
    ## fit the model
    #choose lambda
    #lams <- mcvs01 <- mcvs09 <- seq(0.5, 5, by = 0.5)
    #for(j in 1:length(lams)) {
    #  mcvs01[j] <- MCV(lams[i], y ~ qss(x, lambda = lambdas[1]),
    #                 data.frame(x=x,y=y), tau=0.1, k=5)
    #  mcvs09[j] <- MCV(lams[i], y ~ qss(x, lambda = lambdas[1]),
    #                 data.frame(x=x,y=y), tau=0.9, k=5)
    #}
    #fit the model
    lambda=1.7
    mod1 <- rqss(y ~ qss(x, constraint= "N",lambda = lambda),
                 tau=tau[1])
    
    mod2 <- rqss(y ~ qss(x, constraint= "N",lambda = lambda),
                 tau=tau[2])
    # simulate a new observation to predict
    if(err.distr=="normal"){
      ynew <-beta[1] + beta[2]*(xnew) + beta[3]*(xnew^2) + beta[4]*(xnew^3) + rnorm(1,0,1)
    }
    # compute confidence and prediction intervals
    pi <- c(predict.rqss(mod1, newdata = data.frame(x=xnew)),
            predict.rqss(mod2, newdata = data.frame(x=xnew)))
    if ( (pi[1] < ynew) & (ynew < pi[2])){
      pi.hits <- pi.hits + 1
    }
    pred_bounds[i,] <- pi
    interval_score[i, ] <- interval_score(
      lower_bound = pi[1],
      upper_bound = pi[2],
      observation = ynew,
      alpha = 0.2
    )
  }
  return( list( pi.hits/nsim, pred_bounds, interval_score ) )
}

##simulation coverage 80% quantile regression i.i.d case
set.seed(21)
nibs=1000
nsim=1000
quantreg.coverage.rqss <-  matrix(NA, nsim, 3)
quantreg.intScore.rqss <- matrix(NA, nsim, 3)
cubic_bounds <- list(matrix(NA, nsim, 2), 
                      matrix(NA, nsim, 2),
                      matrix(NA, nsim, 2))

for(k in 1:nsim){
  
  #independent variable
  max.x <- 8
  min.x <- -3
  x <- runif(nobs,min = min.x, max = max.x )
  xnew.values <- qunif(c(.1,.5,.9), min = min.x, max = max.x)
  
  #base model
  beta=c(-0.2, 1, -0.8, 0.12)
  
  for(i in 1:length(xnew.values)){
    
    #normal error
    cubic_pred_int <- coverage.rq.cubic(x=x, 
                                         beta = beta, err.distr = "normal",
                                         tau=c(0.1,0.9), xnew = xnew.values[i],
                                         nsim = 1)
    cubic_bounds[[i]][k,]<- cubic_pred_int[[2]]
    quantreg.coverage.rqss[k,i] <- cubic_pred_int[[1]]
    quantreg.intScore.rqss[k,i] <- cubic_pred_int[[3]]
    
    
  }
  print(k)
}

#lambda=1.7
# coverage
colMeans(quantreg.coverage.rqss)
# [1] 0.774 0.779 0.788

#cubic bounds
mean(cubic_bounds[[1]][,2]-cubic_bounds[[1]][,1])
mean(cubic_bounds[[2]][,2]-cubic_bounds[[2]][,1])
mean(cubic_bounds[[3]][,2]-cubic_bounds[[3]][,1])

#> mean(cubic_bounds[[1]][,2]-cubic_bounds[[1]][,1])
#[1] 2.500651
#> mean(cubic_bounds[[2]][,2]-cubic_bounds[[2]][,1])
#[1] 2.547875
#> mean(cubic_bounds[[3]][,2]-cubic_bounds[[3]][,1])
#[1] 2.515173

# avg interval score
colMeans(quantreg.intScore.rqss)
#  3.562 3.548 3.523


#####

