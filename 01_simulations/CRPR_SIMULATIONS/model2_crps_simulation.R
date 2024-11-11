### CRPS SIMULATION - LINK SELECTION - Model 8 i.i.d student errors(df=3) ###

#####Libraries#####
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

#### CRPS SIMULATION ####

# Settings
nsim <- 250
nobs <- 1000

# Links to be used in VTM 
links <- list("robit1" <- binomial(link = gosset(nu = 1)),
              "robit2" <- binomial(link = gosset(nu = 2)),
              "robit3" <- binomial(link = gosset(nu = 3)),
              "robit4" <- binomial(link = gosset(nu = 4)),
              "robit5" <- binomial(link = gosset(nu = 5)),
              "cloglog" <- binomial(link = "cloglog"),
              "probit" <- binomial(link = "probit"),
              "logit" <- binomial(link = "logit"))

# Matrix to store CRPS results
CRSP_score <- matrix( NA, nsim, length(links) )
colnames(CRSP_score) <- c("robit1",
                          "robit2",
                          "robit3",
                          "robit4",
                          "robit5",
                          "cloglog",
                          "probit",
                          "logit" )

# Generate n samples to train VTM models
independent_var <- matrix( NA, nsim, nobs )
response_var <- matrix( NA, nsim, nobs )
# seed
set.seed(21)
for( k in 1:nsim){

  #independent variable
  mu.x <- 5
  sd.x <- 1
  x <- rnorm(nobs, mean = mu.x, sd =sd.x)
  independent_var[k, ] <- x
  #base model
  beta0=1
  beta1=2
  y=beta0 + beta1*x
  #iid errors chisquare df=3
  df=3
  err_chisq <- rchisq(nobs,df=df)
  #positive chisq model
  y.chisq <- y+err_chisq
  response_var[k, ] <- y.chisq
  
}

# CRPS SIMULATION
# seed
set.seed(21)

for( sim in 1:nsim) {
  for (l in 1:length(links)) {
    
    #dataset
    x <- independent_var[sim, ]
    y <- response_var[sim, ]
    
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
    datp <- data.frame(x = x) ### sample to predict (estimate distr. on train data)
    
    # VTM
    Splitfit <- VTM(
      data_train = ds,
      data_predict = datp,
      formula = formbin,
      splits = splits,
      num_interpolation_points = length(splits) * 2,
      family_link = links[[l]],
      lambda = 0,
      alpha = 0
    )
    
    # CRPS
    mean_crps <- c()
    for(i in 1:length(y)){
      mean_crps <- c( mean_crps,  
                      crps_sample(y[i],
                                  dat=Splitfit$interpolation_points,
                                  w = Splitfit$histogram[i,]) )
    }
    CRSP_score[sim, l] <- mean(mean_crps)
  }
  print(sim)
}

table( colnames(CRSP_score)[ apply(CRSP_score[1:nsim,], 1, FUN = which.min) ] )
#cloglog 
#    250

table( colnames(CRSP_score)[ apply(CRSP_score[1:nsim,], 1, FUN = which.min) ] )/ nsim
#cloglog 
#    1


write.csv(x=CRSP_score,file = "C:\\Users\\nicco\\Desktop\\UNI\\00Esami in corso\\Tesi\\ADAC_revision\\crps_score.csv")
