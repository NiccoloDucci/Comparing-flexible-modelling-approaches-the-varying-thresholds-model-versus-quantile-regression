##### functions#####

#### START FUNCTION
VTM <- function(data_train,
                data_predict,
                formula,
                splits,
                num_interpolation_points,
                family_link,
                lambda,
                alpha) {
  ## input
  # data_train:  learning data
  # datap:  prediction data
  # splits: splits in [min, max] of response variable range
  # formula: regression formula (in formula response has name binresp, e.g., binresp ~ x )
  # num_interpolation_points:  number of x-values to compute distribution_function and histogram.
  # family_link: link function, e.g., "binomial(link = "logit")", if alpha or lambda > 0 family_link = binomial()
  # lambda: >0 ridge
  # alpha: (lasso paramter), elasticnet possible,  alpha =1 (lasso) possible, alpha=0 (ridge)
  ##

  ## minimum observations
  minobsfit <- 4

  ## link function
  family1 <- family_link

  ## regression settings
  funct_form <- gsub(" ", "", formula[3])
  groupvars <- unlist(strsplit(funct_form, "\\+"))
  numvar <- length(groupvars) + 1
  numsplits <- length(splits)
  n_rows_data_train <- dim(data_train)[1]
  param <- matrix(0, numvar, numsplits)
  param1 <- matrix(0, numvar, 1)
  stderr <- matrix(0, numvar, numsplits)
  n_rows_data_predict <- dim(data_predict)[1]

  ## response values for which distribution probability will be computed
  minobs <- min(data_train$resp)
  maxobs <- max(data_train$resp)
  interpolation_points <- seq(minobs - 0.001, maxobs, (maxobs - minobs) / num_interpolation_points)
  legnth_inter_points <- length(interpolation_points)

  ## vectors where to store estimated probabilities
  distribution_function_splitsdum <- matrix(0, n_rows_data_predict, numsplits)
  distribution_function_splits <- matrix(0, n_rows_data_predict, numsplits)

  ###### loop for distribution_function on split extremes
  valid <- 0

  ################ glm fit
  if (lambda <= 0) {
    for (l in 1:numsplits) {
      data_train$binresp <- 0
      for (i in 1:n_rows_data_train) {
        if (data_train$resp[i] > splits[l]) {
          data_train$binresp[i] <- 1
        }
      }

      ## number of observation for 0s and 1s
      nobs0 <- sum(data_train$binresp == 0)
      nobs1 <- sum(data_train$binresp == 1)

      ###### if for fit
      if ((nobs0 > minobsfit) & (nobs1 > minobsfit)) {
        valid <- c(valid, l)
        glmfit <- glm(formula, data = data_train, family = family1) # qui21
        param[, l] <- glmfit$coefficients
        stderr[, l] <- matrix(sqrt(diag(vcov(glmfit))), numvar, 1)
        d <- predict(glmfit, data_predict, type = "response")
        dumd <- d
        distribution_function_splitsdum[, l] <- 1 - d
      } # end if for fit
    }

    ### fit extremes
    valid <- valid[2:length(valid)]
    minv <- min(valid)
    maxv <- max(valid)
    if (minv > 1) {
      data_train$binresp <- 0
      for (i in 1:n_rows_data_train) {
        if (data_train$resp[i] > splits[minv]) {
          data_train$binresp[i] <- 1
        }
      }
      glmfit <- glm(formula, data = data_train, family = family1) # qui21
      minv1 <- minv - 1
      for (li in 1:minv1) {
        param[, li] <- glmfit$coefficients
        stderr[, li] <- matrix(sqrt(diag(vcov(glmfit))), numvar, 1)
        d <- predict(glmfit, data_predict, type = "response")
        dumd <- d
        distribution_function_splitsdum[, li] <- 1 - d
      }
    } # end if minv
    if (maxv < numsplits) {
      data_train$binresp <- 0
      for (i in 1:n_rows_data_train) {
        if (data_train$resp[i] > splits[maxv]) {
          data_train$binresp[i] <- 1
        }
      }
      glmfit <- glm(formula, data = data_train, family = family1) # qui21
      start <- maxv + 1
      for (li in start:numsplits) {
        param[, li] <- glmfit$coefficients
        stderr[, li] <- matrix(sqrt(diag(vcov(glmfit))), numvar, 1)
        d <- predict(glmfit, data_predict, type = "response")
        dumd <- d
        distribution_function_splitsdum[, li] <- 1 - d
      }
    } # end if max
  }
  ### end if lambda <=0

  ### ridge/lasso fit
  if (lambda > 0 | alpha > 0) {
    for (l in 1:numsplits) {
      data_train$binresp <- 0
      for (i in 1:n_rows_data_train) {
        if (data_train$resp[i] > splits[l]) {
          data_train$binresp[i] <- 1
        }
      }
      data_train$binresp <- as.factor(data_train$binresp)

      ##### glmnet
      x_var <- data.matrix(data_train[, groupvars])
      y_var <- data_train[, "binresp"]
      netfit <- glmnet(x_var, y_var, family = binomial(), alpha = alpha, lambda = lambda)
      param1 <- as.numeric(coef(netfit))
      param[, l] <- param1
      ######

      ## prediction new
      newx <- data.matrix(data_predict[, groupvars])
      pr <- predict(netfit, newx, s = lambda, type = c("response"))
      distribution_function_splitsdum[, l] <- 1 - pr
    }
  }
  ### end if ridge/lasso

  ## monotone regression to fix parameters
  for (i in 1:n_rows_data_predict) {
    iso <- isoreg(splits, distribution_function_splitsdum[i, ])
    distribution_function_splits[i, ] <- iso$yf
  }

  ## vectors for estimated probabilities
  distribution_function_jump <- matrix(0, nrow = n_rows_data_predict, ncol = legnth_inter_points)
  distribution_function_poly <- matrix(0, nrow = n_rows_data_predict, ncol = legnth_inter_points)

  ## distribution function no interpolation
  for (i in 1:n_rows_data_predict) {
    for (j in 1:numsplits) {
      for (l in 1:legnth_inter_points) {
        if (interpolation_points[l] > splits[j]) {
          distribution_function_jump[i, l] <- distribution_function_splits[i, j]
        }
      }
    }
  }

  ## distribution function with interpolation
  ones <- matrix(1, n_rows_data_predict, 1)
  distribution_function_splitsapp <- cbind(distribution_function_splits, ones)
  distsplits <- splits[2] - splits[1]

  for (i in 1:n_rows_data_predict) {
    for (l in 1:legnth_inter_points) {
      if (interpolation_points[l] > minobs) {
        distribution_function_poly[i, l] <-
          (distribution_function_splitsapp[i, 1] / (splits[1] - minobs)) * (interpolation_points[l] - minobs)
      }
    }
    numsplits1 <- numsplits - 1
    for (j in 1:numsplits1) {
      for (l in 1:legnth_inter_points) {
        if (interpolation_points[l] > splits[j]) {
          distribution_function_poly[i, l] <-
            distribution_function_splitsapp[i, j] + ((distribution_function_splitsapp[i, j + 1] - distribution_function_splitsapp[i, j]) / distsplits) * (interpolation_points[l] - splits[j])
        }
      }
    }
    for (l in 1:legnth_inter_points) {
      if (interpolation_points[l] > splits[numsplits]) {
        distribution_function_poly[i, l] <-
          distribution_function_splitsapp[i, numsplits] +
          ((1 - distribution_function_splitsapp[i, numsplits]) / (maxobs - splits[numsplits])) * (interpolation_points[l] - splits[numsplits])
      }
    }

    for (l in 1:legnth_inter_points) {
      if (interpolation_points[l] > maxobs) distribution_function_poly[i, l] <- 1
    }
  }

  ## final distribution function with interpolation
  distribution_function <- distribution_function_poly



  #### histogram (pdf) - check negative values due to small (1e-10) rounding errors
  histp <- matrix(0, n_rows_data_predict, legnth_inter_points)
  for (i in 1:n_rows_data_predict) {
    for (l in 2:legnth_inter_points) {
      pdf_diff_estimate <- distribution_function_poly[i, l] - distribution_function_poly[i, l - 1]
      if (pdf_diff_estimate >= 0) {
        histp[i, l] <- pdf_diff_estimate
      } else if (abs(pdf_diff_estimate) < 1e-10) {
        histp[i, l] <- 0
      } else {
        warning("negative values in pdf estimation")
      }
    }
  }

  ### mean
  meanval <- matrix(0, n_rows_data_predict, 1)
  for (i in 1:n_rows_data_predict) {
    for (l in 1:legnth_inter_points) meanval[i, 1] <- meanval[i, 1] + interpolation_points[l] * histp[i, l]
  }

  ### median
  medianval <- matrix(0, n_rows_data_predict, 1)
  numxax1 <- legnth_inter_points - 1
  for (i in 1:n_rows_data_predict) {
    for (l in 1:numxax1) if (distribution_function[i, l] <= 0.5) medianval[i, 1] <- (interpolation_points[l] + interpolation_points[l + 1]) / 2
  }


  newList <- list(
    "histogram" = histp,
    "distribution_function" = distribution_function,
    "parameters" = param,
    "stderr" = stderr,
    "interpolation_points" = interpolation_points,
    "distribution_function_jump" = distribution_function_jump,
    "meanval" = meanval,
    "medianval" = medianval
  )

  ## output
  # histogram: estimated probability density function
  # distribution_function: estimated distribution function
  # parameters: binary models parameters
  # stderr: standard deviation of binary models parameters
  # interpolation_points: response variable values used for interpolation between thresholds
  # distribution_function_jump: estimated distribution function using only prior thresholds
  # meanval: vector of conditional means given the covariate value
  # medianval: vector of conditional medians given the covariate value
  ##

  return(newList)
}


################## end function
