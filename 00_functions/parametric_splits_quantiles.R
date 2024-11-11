
#### quantile function

ParametricSplitsQuantiles <- function(distribution_function, y_values, quantile) {
  ## input
  # distribution_function: matrix with distribution_function
  # quantile: quantile
  # y_values: values where distribution_function is computed (ideally equally distribuited within the support)
  ##
  
  n_rows <- dim(distribution_function)[1]
  n_columns <- dim(distribution_function)[2]
  quantilevalues <- matrix(0, n_rows, 1)
  
  for (i in 1:n_rows) {
    for (j in 1:n_columns) {
      if (distribution_function[i, j] <= quantile) {
        quantilevalues[i, 1] <- y_values[j]
      }
    }
    quantilevalues[i, 1] <- quantilevalues[i, 1] + (y_values[2] - y_values[1]) / 2
  }
  
  ## output
  # quantilevalues: quantile values estimated
  ##
  
  newList <- list("quantilevalues" = quantilevalues)
  return(newList)
}

#### end quantile function
