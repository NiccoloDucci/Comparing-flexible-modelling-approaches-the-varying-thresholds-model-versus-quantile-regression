### PREDICTION INTERVALS - model 8 t student (df=3) i.i.d. errors with quadratic effect ###

##### Libraries#####
library(tidyverse)
library(ggplot2)
library(quantreg)
#####

##### functions #####

source("C:\\Users\\nicco\\Desktop\\UNI\\00Esami in corso\\Tesi\\ADAC_revision\\00_functions\\varying_thresholds_model.R")
source("C:\\Users\\nicco\\Desktop\\UNI\\00Esami in corso\\Tesi\\ADAC_revision\\00_functions\\parametric_splits_quantiles.R")
source("C:\\Users\\nicco\\Desktop\\UNI\\00Esami in corso\\Tesi\\ADAC_revision\\00_functions\\interval_score.R")

#####


#### model 1 normal covariate hairplot ####

### SETTINGS
nsim <- 70
nobs <- 1000

# seed
set.seed(211)

# arrays for quantiles
Splitfit.quant10.norm <- matrix(NA, nrow = 1000, ncol = nsim)
Splitfit.quant90.norm <- matrix(NA, nrow = 1000, ncol = nsim)



for (i in 1:nsim) {
  # independent variable
  mu.x <- 5
  sd.x <- 1
  x <- rnorm(nobs, mean = mu.x, sd = sd.x)

  # iid errors standard normal
  err_normal <- rnorm(nobs, 0, 1)

  # base model
  beta0 <- 1
  beta1 <- 2
  y <- beta0 + beta1 * x

  # normal model
  y.response <- y + err_normal

  ## dataset
  datset <- data.frame(y = y.response, x = x)


  ## varying-threshold model

  ### fit the model
  # splits generation
  # k=52 (50 + 2), there are k+1 theta thresholds
  splitnumber <- 50
  # theta1 - thetak-1
  minsplit <- quantile(datset$y, probs = 0.02)
  maxsplit <- quantile(datset$y, probs = 0.98)
  increase <- (maxsplit - minsplit) / splitnumber
  splits <- seq(minsplit, maxsplit, by = increase)
  # formula
  formbin <- binresp ~ x
  # dataset modification
  ds <- data.frame(y = datset$y, x = datset$x)
  names(ds)[1] <- "resp"
  datp <- data.frame(x = datset$x) ### sample to predict

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



  # estimated VTM quantiles
  Splitfit.quant10.norm[, i] <- ParametricSplitsQuantiles(
    distribution_function = Splitfit$distribution_function,
    quantile = .1, y_values = Splitfit$interpolation_points
  )$quantilevalues %>% sort()

  Splitfit.quant90.norm[, i] <- ParametricSplitsQuantiles(
    distribution_function = Splitfit$distribution_function,
    quantile = .9, y_values = Splitfit$interpolation_points
  )$quantilevalues %>% sort()

  print(i)
}


# x aixs
x.axis <- qnorm(seq(0.0001, 0.9999, 0.001), 0, 10)

# data frame of quantile lines
joint_quant_df <- cbind(Splitfit.quant10.norm, Splitfit.quant90.norm)

# matplot
matplot(
  x = x.axis, y = joint_quant_df[, c(1:70, 71:140)],
  type = "n",
  ylab = "y",
  xlab = "x",
  xaxt = "n",
  yaxt = "n",
  cex.lab = 2.4
) # this first matplot allows the grid to not overlaps the next matplot
grid(col = "grey90",
     lty = 1
)
matplot(
  x = x.axis, y = joint_quant_df[, c(1:70, 71:140)], type = "l", col = "black",
  lty = 1, ylab = "y", xlab = "x", add = T
)

# Customize text colors
title(col.main = "grey20", col.lab = "grey20") # Grey title and axis labels
axis(1, col.axis = "grey20", cex.axis = 2) # Grey x-axis text
axis(2, col.axis = "grey20", cex.axis = 2) # Grey y-axis text

#####


#### model 1 uniform covariate hairplot ####

### SETTINGS
nsim <- 70
nobs <- 1000

# seed
set.seed(211)

# arrays for quantiles
Splitfit.quant10.unif <- matrix(NA, nrow = 1000, ncol = nsim)
Splitfit.quant90.unif <- matrix(NA, nrow = 1000, ncol = nsim)



for (i in 1:nsim) {
  
  # independent variable
  x <- runif(nobs,2, 8)
  
  # iid errors standard normal
  err_normal <- rnorm(nobs, 0, 1)
  
  # base model
  beta0 <- 1
  beta1 <- 2
  y <- beta0 + beta1 * x
  
  # normal model
  y.response <- y + err_normal
  
  ## dataset
  datset <- data.frame(y = y.response, x = x)
  
  
  ## varying-threshold model
  
  ### fit the model
  # splits generation
  # k=52 (50 + 2), there are k+1 theta thresholds
  splitnumber <- 50
  # theta1 - thetak-1
  minsplit <- quantile(datset$y, probs = 0.02)
  maxsplit <- quantile(datset$y, probs = 0.98)
  increase <- (maxsplit - minsplit) / splitnumber
  splits <- seq(minsplit, maxsplit, by = increase)
  # formula
  formbin <- binresp ~ x
  # dataset modification
  ds <- data.frame(y = datset$y, x = datset$x)
  names(ds)[1] <- "resp"
  datp <- data.frame(x = datset$x) ### sample to predict
  
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
  
  
  
  # estimated VTM quantiles
  Splitfit.quant10.unif[, i] <- ParametricSplitsQuantiles(
    distribution_function = Splitfit$distribution_function,
    quantile = .1, y_values = Splitfit$interpolation_points
  )$quantilevalues %>% sort()
  
  Splitfit.quant90.unif[, i] <- ParametricSplitsQuantiles(
    distribution_function = Splitfit$distribution_function,
    quantile = .9, y_values = Splitfit$interpolation_points
  )$quantilevalues %>% sort()
  
  print(i)
}


# x aixs
x.axis <- qunif(seq(0.0001, 0.9999, 0.001), 2, 8)

# data frame of quantile lines
joint_quant_df.unif <- cbind(Splitfit.quant10.unif, Splitfit.quant90.unif)

# matplot
matplot(
  x = x.axis, y = joint_quant_df.unif[, c(1:70, 71:140)],
  type = "n",
  ylab = "y",
  xlab = "x",
  xaxt = "n",
  yaxt = "n",
  cex.lab = 2.4
) # this first matplot allows the grid to not overlaps the next matplot
grid(col = "grey90",
     lty = 1
)
matplot(
  x = x.axis, y = joint_quant_df.unif[, c(1:70, 71:140)], type = "l", col = "black",
  lty = 1, ylab = "y", xlab = "x", add = T
)

# Customize text colors
title(col.main = "grey20", col.lab = "grey20") # Grey title and axis labels
axis(1, col.axis = "grey20", cex.axis = 2,col.ticks = "grey20") # Grey x-axis text
axis(2, col.axis = "grey20", cex.axis = 2,col.ticks = "grey20") # Grey y-axis text

#####


#### combined charts ####

par(mfrow = c(1, 2))

### normal chart
# x aixs
x.axis <- qnorm(seq(0.0001, 0.9999, 0.001), 5, 1)

# data frame of quantile lines
joint_quant_df <- cbind(Splitfit.quant10.norm, Splitfit.quant90.norm)

# matplot
matplot(
  x = x.axis, y = joint_quant_df[, c(1:70, 71:140)],
  type = "n",
  ylab = "y",
  xlab = "x",
  xaxt = "n",
  yaxt = "n",
  cex.lab = 2.4,
  ylim = c(0, 20)
) # this first matplot allows the grid to not overlaps the next matplot
grid(col = "grey90",
     lty = 1
)
matplot(
  x = x.axis, y = joint_quant_df[, c(1:70, 71:140)], type = "l", col = "black",
  lty = 1, ylab = "y", xlab = "x", add = T
)

# Customize text colors
title(col.main = "grey20", col.lab = "grey20") # Grey title and axis labels
axis(1, col.axis = "grey20", cex.axis = 2, tck=-0.01) # Grey x-axis text
axis(2, col.axis = "grey20", cex.axis = 2, tck=-0.01) # Grey y-axis text
###




### uniform chart
# x aixs
x.axis <- qunif(seq(0.0001, 0.9999, 0.001), 2, 8)

# data frame of quantile lines
joint_quant_df.unif <- cbind(Splitfit.quant10.unif, Splitfit.quant90.unif)

# matplot
matplot(
  x = x.axis, y = joint_quant_df.unif[, c(1:70, 71:140)],
  type = "n",
  ylab = "y",
  xlab = "x",
  xaxt = "n",
  yaxt = "n",
  cex.lab = 2.4,
  ylim = c(0, 20)
) # this first matplot allows the grid to not overlaps the next matplot
grid(col = "grey90",
     lty = 1
)
matplot(
  x = x.axis, y = joint_quant_df.unif[, c(1:70, 71:140)], type = "l", col = "black",
  lty = 1, ylab = "y", xlab = "x", add = T
)

# Customize text colors
title(col.main = "grey20", col.lab = "grey20") # Grey title and axis labels
axis(1, col.axis = "grey20", cex.axis = 2, tck=-0.01) # Grey x-axis text
axis(2, col.axis = "grey20", cex.axis = 2, tck=-0.01) # Grey y-axis text
###

#####