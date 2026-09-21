###############################################
# ramey-zubairy-2018-lp
# Author- Amal Varghese
#  Purpose: LP function
###############################################


# Prelim ------------------------------------------------------------------

rm(list=ls())
library(dplyr); library(sandwich)

rzdat<- readRDS("data/processed/rzdat_prepped.rds")
glimpse(rzdat)


# Generate lagged controls ---------------------------------------------------

make_lag_matrix<- function(data, vars, p) {
  var_blocks <- lapply(vars, function(variable){
    block <- sapply (1:p, function (lag_number) dplyr::lag(data[[variable]], lag_number))
    colnames(block) <- paste0(variable, "_L", 1:p)
    block
  })
  do.call(cbind, var_blocks)
}

lag_controls <- make_lag_matrix(rzdat, vars = c("newsy", "y", "g"), p = 4)
 dim(lag_controls)
 

# The single-horizon regression,  ℎ = 0 -----------------------------------

h <- 0
y_h <- dplyr::lead(rzdat$y, h)

# GDP response 
fit0 <- lm(y_h ~ rzdat$newsy  + lag_controls)
summary(fit0)

# government spending response
g_h <- dplyr::lead(rzdat$g, h)

fit0_g <- lm(g_h ~ rzdat$newsy + lag_controls)
summary(fit0_g)


# NW Standard Errors ------------------------------------------------------

library(lmtest)

nw_vcov_y <- NeweyWest(fit0, lag = 4, prewhite = FALSE, adjust = TRUE)
coeftest(fit0, vcov = nw_vcov_y)

nw_vcov_g <- NeweyWest(fit0_g, lag = 4, prewhite = FALSE, adjust = TRUE)
coeftest(fit0_g, vcov = nw_vcov_g)




# LP FUNCTION -------------------------------------------------------------

estimate_horizon <- function(h, outcome, shock, controls, data, hac_lag = 4) {
  y_h <- dplyr::lead(data[[outcome]], h)
  fit <- lm(y_h ~ data[[shock]] + controls)
  
  nw_vcov <- NeweyWest(fit, lag = hac_lag, prewhite = FALSE, adjust = TRUE)
  ct <- coeftest(fit, vcov = nw_vcov)
  
  shock_row <- ct["data[[shock]]", ]
  
  data.frame(
    h = h,
    coef = shock_row["Estimate"],
    se = shock_row["Std. Error"],
    row.names = NULL
  )
}

# trial on h=0
estimate_horizon(h = 0, outcome = "y", shock = "newsy", controls = lag_controls, data = rzdat)
estimate_horizon(h = 0, outcome = "g", shock = "newsy", controls = lag_controls, data = rzdat)


# Looping for IRF ---------------------------------------------------------

library(purrr)

irf_y_linear <- map_dfr(0:20, ~ estimate_horizon(
  h = .x, outcome = "y", shock = "newsy", controls = lag_controls, data = rzdat
))

irf_g_linear <- map_dfr(0:20, ~ estimate_horizon(
  h = .x, outcome = "g", shock = "newsy", controls = lag_controls, data = rzdat
))

irf_y_linear
irf_g_linear



# Plot Linear IRFs --------------------------------------------------------

library(ggplot2)
library(patchwork)

plot_irf <- function(irf_df) {
  ggplot(irf_df, aes(x = h, y = coef)) +
    geom_ribbon(aes(ymin = coef - 1.96 * se, ymax = coef + 1.96 * se), alpha = 0.2) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      aspect.ratio = 1
    )
}

p_g <- plot_irf(irf_g_linear)
p_y <- plot_irf(irf_y_linear)

p_g / p_y

dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)
ggsave("output/figures/figure5_linear.png", p_g / p_y, width = 6, height = 8, dpi = 150)