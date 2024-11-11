### Interval score from Gneiting and Raftery (2007)

interval_score <- function(lower_bound,
                           upper_bound,
                           observation,
                           alpha) {
  interval <- (upper_bound - lower_bound)
  penalty <- (2 / alpha * (lower_bound - observation)) * (observation < lower_bound) * 1 + (2 / alpha * (observation - upper_bound)) * (observation > upper_bound) * 1
  score <- interval + penalty
  return(score)
}

###
