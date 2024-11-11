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


#####


#### SIMULATED DATASET MODEL 5 ####

### SETTINGS
nobs <- 1000

# seed
set.seed(21)

# independent variable
x <- runif(nobs,-2,12)

# iid errors standard normal
err_normal <- rnorm(nobs, 0, 1)

# base model
beta0=0.2
beta1=-0.4
beta2=0.1
y=beta0 + beta1*x +beta2*x^2 

# model 5
y.response <- y + err_normal

## dataset
datset <- data.frame(y = y.response, x = x)

#####


#### quantile regression ####

## quantile regression
taus <- seq(.01, .99, .01)
quantreg <- rq(formula = y ~ x, data = datset, tau = taus)

#####

#### varying-threshold model ####

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

#####


#### scatter plot ####

# estimated VTM quantiles
Splitfit.quant10 <- ParametricSplitsQuantiles(
  distribution_function = Splitfit$distribution_function,
  quantile = .1, y_values = Splitfit$interpolation_points
)

Splitfit.quant50 <- ParametricSplitsQuantiles(
  distribution_function = Splitfit$distribution_function,
  quantile = .5, y_values = Splitfit$interpolation_points
)

Splitfit.quant90 <- ParametricSplitsQuantiles(
  distribution_function = Splitfit$distribution_function,
  quantile = .9, y_values = Splitfit$interpolation_points
)

# x and y aixs ticks labels
min.x <- min(datset$x)
max.x <- max(datset$x)
sd.x <- sd(datset$x)
min.y <- min(datset$y)
max.y <- max(datset$y)
sd.y <- sd(datset$y)
x_axis <- seq(floor(min.x), ceiling(max.x), round(sd.x, 0))
y_axis <- seq(floor(min.y), ceiling(max.y), round(sd.y, 0))



# scatter plot
ggplot() +
  ggtitle("") +
  labs(color = "Legend:") +
  geom_vline(
    xintercept = x_axis,
    linetype = 1, color = "grey90"
  ) +
  geom_hline(
    yintercept = y_axis,
    linetype = 1, color = "grey90"
  ) +
  geom_point(aes(x = datset$x, y = datset$y),
             shape = 18,
             size = 2
  ) +
  geom_line(
    aes(
      x = datset$x,
      y = quantreg$fitted.values[, 10],
      color = "Quantile regression"
    ),
    linewidth = 1.5
  ) +
  geom_line(
    aes(
      x = datset$x,
      y = quantreg$fitted.values[, 50]
    ),
    linewidth = 1.5, color = "grey60"
  ) +
  geom_line(
    aes(
      x = datset$x,
      y = quantreg$fitted.values[, 90]
    ),
    linewidth = 1.5, color = "grey60"
  ) +
  geom_line(
    aes(
      x = datset$x,
      y = Splitfit.quant10$quantilevalues,
      color = "VTM"
    ),
    linewidth = 1.5
  ) +
  geom_line(
    aes(
      x = datset$x,
      y = Splitfit.quant50$quantilevalues
    ),
    linewidth = 1.5, color = "tomato2"
  ) +
  geom_line(
    aes(
      x = datset$x,
      y = Splitfit.quant90$quantilevalues
    ),
    linewidth = 1.5, color = "tomato2"
  ) +
  scale_color_manual(
    values = c("VTM" = "tomato2", "Quantile regression" = "grey60"),
    labels = c("Quantile regression", "VTM") # Custom labels for the legend
  ) +
  theme_classic() +
  scale_x_continuous(
    name = "x",
    breaks = x_axis
  ) +
  scale_y_continuous(
    name = "y",
    breaks = y_axis
  ) +
  theme(
    axis.text.x = element_text(
      size = 25,
      colour = "grey20"
    ),
    axis.text.y = element_text(
      size = 25,
      colour = "grey20"
    ),
    axis.title = element_text(
      size = 30,
      face = "plain",
      colour = "grey20"
    ),
    plot.title = element_text(
      hjust = 0.5,
      size = 20
    ),
    legend.position = "top",
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 18),
    panel.border = element_rect(color = "black", fill = NA, size = 1), # Add black border around plot area
    plot.margin = margin(5, 10, 5, 5), # Optional: Adjust margins to avoid clipping)
    legend.key.width = unit(2, "cm"), # Increase width of legend key (line length)
    legend.key.height = unit(1, "cm")
  )


#####


#### VTM binary regressions parameters plot ####

# settings
param <- Splitfit$parameters
stderr <- Splitfit$stderr
par_number <- seq(0, nrow(param), 1) %>% as.character()
par_color <- c("black", grey.colors(nrow(param) - 1))
matplot_lables <- paste0("split", "_", 1:length(splits), ": ", round(splits[1:length(splits)], 1))

# all variables
label.size <- 1.6
matplot(splits, t(param),
        type = "b",
        main = "VTM Parameter Functions: coefficients vs splits",
        xlab = "",
        ylab = "Varying beta parameters",
        xaxt = "n",
        lwd = 3,
        cex.lab = label.size,
        cex.axis = label.size,
        cex = 1.6,
        cex.main = 1.5,
        pch = par_number,
        col = par_color
)
abline(v = splits, lty = "dashed", col = "tomato2")
abline(h = 0, lty = "dashed", col = "black", lwd = 1.3)
axis(1, at = splits, labels = FALSE)
text(x = splits, y = par("usr")[3] - 0.4, labels = matplot_lables, srt = 45, xpd = TRUE, adj = 1, cex = 1.15)
#####



#### comparison between varying beta ####

# settings
param <- Splitfit$parameters
stderr <- Splitfit$stderr
var <- 2 # which variable to plot(intercept or covariates)
taus <- seq(.01, .99, .01)
quantiles_y <- quantile(datset$y, taus)
# x and y aixs ticks labels
min.y <- min(quantiles_y)
max.y <- max(quantiles_y)
sd.y <- sd(quantiles_y)
y_axis <- seq(floor(min.y), ceiling(max.y), round(sd.y, 0))

comparison_plot <- ggplot() +
  ggtitle("") +
  ylab("Regression Coefficients") +
  xlab("Response variable") +
  labs(color = "Legend:") +
  geom_vline(
    xintercept = y_axis,
    linetype = 1, color = "grey90"
  ) +
  geom_hline(
    yintercept = seq(
      0,
      ceiling(max(
        param[var, ] + 1.96 * stderr[var, ],
        floor(quantreg$coefficients[var, ])
      )),
      1
    ),
    linetype = 1, color = "grey90"
  ) +
  geom_step(
    aes(
      x = splits,
      y = param[2, ],
      color = "VTM"
    ),
    linewidth = 1.5
  ) +
  geom_step(
    aes(
      x = quantiles_y,
      y = quantreg$coefficients[2, ],
      color = "Quantile regression"
    ),
    linewidth = 1.5
  ) +
  geom_hline(aes(yintercept = beta1, col = "True beta"),
             linewidth = 1.3,
             linetype = 2
  ) +
  scale_color_manual(
    values = c("VTM" = "tomato2", "Quantile regression" = "grey60", "True beta" = "black"),
    labels = c("Quantile regression", "True beta", "VTM") # Custom labels for the legend
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      size = 20,
      colour = "grey20"
    ),
    axis.text.y = element_text(
      size = 20,
      colour = "grey20"
    ),
    axis.title = element_text(
      size = 20,
      colour = "grey20"
    ),
    plot.title = element_text(hjust = 0.5, size = 20),
    legend.position = "top",
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 18),
    panel.border = element_rect(color = "black", fill = NA, size = 1), # Add black border around plot area
    plot.margin = margin(5, 10, 5, 5), # Optional: Adjust margins to avoid clipping)
    legend.key.width = unit(2, "cm"), # Increase width of legend key (line length)
    legend.key.height = unit(1, "cm")
  )

comparison_plot +
  geom_line(aes(x = splits, y = param[var, ] + 1.96 * stderr[var, ]),
            linetype = 2,
            linewidth = 0.8,
            col = "tomato2"
  ) +
  geom_line(aes(x = splits, y = param[var, ] - 1.96 * stderr[var, ]),
            linetype = 2,
            linewidth = 0.8,
            col = "tomato2"
  ) +
  scale_x_continuous(breaks = y_axis)


#####
