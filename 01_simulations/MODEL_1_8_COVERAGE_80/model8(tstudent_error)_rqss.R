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
set.seed(202)

#independent variable
nobs=1000
x <- runif(nobs,0,6)

#iid errors student df=3
df=3
err_student <- rt(nobs,df=df)

#base model
beta0=1
beta1=0.5
beta2=1
y=beta0 + beta1*x + beta2*(x^2)

#positive student model
y.student <- y+err_student


#student
datset.student <- data.frame(y=y.student, x=x)

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
                 datset.student)
par(mfrow = c(1,2))
plot(lams, mcvs, cex = .5, lwd = 2, type = 'l',
     xlab = expression(lambda), ylab = expression(MCV( lambda )))
lambdastar <- lams[which.min(mcvs)]

plot(datset.student$x, datset.student$y, cex = .5, col = "grey")
f <- rqss(y ~ qss(x, lambda = lambdastar),
          data = datset.student, tau=0.5)
plot(f, add = TRUE, lwd = 2)
text(10, 1,bquote(lambda == ~ .(lambdastar)))

#lambda
lambda=lambdastar

#rq splines
mod1 <- rqss(y ~ qss(x, constraint= "N",lambda = lambda), 
             data = datset.student,
             tau=0.1)

mod2 <- rqss(y ~ qss(x, constraint= "N",lambda = lambda), 
             data = datset.student,
             tau=0.9)

ggplot()+
  geom_point(aes(x=datset.student$x, y=datset.student$y))+
  geom_line(aes(x=datset.student$x, 
                y= predict.rqss(mod1,datset.student)),
            col="red",
            lwd=1.6)+
  geom_line(aes(x=datset.student$x, 
                y=predict.rqss(mod2,datset.student)),
            col="red",
            lwd=1.6)+
  theme_classic()
#####


#####GOSSET PREDICTION INTERVAL SIMULATION - student t #####

##quantile regression
coverage.rq <- function(
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
    if(err.distr=="student"){
      err <- rt(n, df=3)
    }
    # simulate the observed data
    y <- beta[1] + beta[2]*x + beta[3]*(x^2) + err
    #fit the model
    lambda=0.7
    mod1 <- rqss(y ~ qss(x, constraint= "N",lambda = lambda),
                 tau=tau[1])
    
    mod2 <- rqss(y ~ qss(x, constraint= "N",lambda = lambda),
                 tau=tau[2])
    # simulate a new observation to predict
    if(err.distr=="student"){
      ynew <- beta[1] + beta[2]*xnew + beta[3]*(xnew^2) + rt(1, df=3)
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
nsim <- 1000
nobs <- 1000
quantreg.coverage.student <- matrix(NA, nsim, 3)
splitfit.intScore.student <- matrix(NA, nsim, 3)
student_bounds <-  list(matrix(NA, nsim, 2), 
                      matrix(NA, nsim, 2),
                      matrix(NA, nsim, 2))

for(k in 1:nsim){
  
  #independent variable
  x <- runif(nobs,0,6)
  xnew.values <- qunif(c(.1,.5,.9), 0,6)
  
  #base model
  beta0=1
  beta1=0.5
  beta2=1
  
  for(i in 1:length(xnew.values)){
    
    #student error
    student_pred_int <- coverage.rq(x=x,
                                  beta = c(beta0,beta1,beta2), err.distr = "student",
                                  tau=c(0.1,0.9), xnew = xnew.values[i],
                                  nsim = 1)
    student_bounds[[i]][k,]<- student_pred_int[[2]]
    quantreg.coverage.student[k,i] <- student_pred_int[[1]]
    splitfit.intScore.student[k,i] <- student_pred_int[[3]]
  }
  print(k)
}

#lambda=0.7
# coverage
colMeans(quantreg.coverage.student)
#[1] 0.770 0.800 0.797

#student bounds
mean(student_bounds[[1]][,2]-student_bounds[[1]][,1])
mean(student_bounds[[2]][,2]-student_bounds[[2]][,1])
mean(student_bounds[[3]][,2]-student_bounds[[3]][,1])
#> mean(student_bounds[[1]][,2]-student_bounds[[1]][,1])
#[1] 3.279746
#> mean(student_bounds[[2]][,2]-student_bounds[[2]][,1])
#[1] 3.267111
#> mean(student_bounds[[3]][,2]-student_bounds[[3]][,1])
#[1] 3.278878

# avg interval score
colMeans(splitfit.intScore.student)
#[1] 5.819894 5.629576 5.557967


#####

