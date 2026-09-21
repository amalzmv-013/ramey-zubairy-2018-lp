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

#save
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)
ggsave("output/figures/figure5_linear.png", p_g / p_y, width = 6, height = 8, dpi = 150)



# Build the state-interacted regressors -----------------------------------

state_lag <- dplyr::lag(rzdat$slack, 1)

shock_high <- rzdat$newsy * state_lag
shock_low <- rzdat$newsy * (1 - state_lag)

controls_high <- lag_controls * state_lag
controls_low <- lag_controls * (1- state_lag)

length(state_lag)
dim(controls_high); dim(controls_low)
table(state_lag, useNA = "always")



#The h=0 state-dependent regression --------------------------------------

h <- 0
y_h <- dplyr::lead(rzdat$y, h)

fit0_state <- lm(y_h ~ 0 + state_lag +
                   shock_high + shock_low + controls_high + controls_low)
summary(fit0_state)

# NW Std Errors

nw_vcov_state <- NeweyWest(fit0_state, lag = 4, prewhite = FALSE, adjust = TRUE)
ct_state <- coeftest(fit0_state, vcov = nw_vcov_state)

ct_state["shock_high", ]
ct_state["shock_low", ]


# SDLP Function -----------------------------------------------------------

estimate_statelp_horizon <- function(h, outcome, state_dummy, shock_high, shock_low, 
                                     controls_high, controls_low, data, hac_lag = 4) {
  y_h <- dplyr::lead(data[[outcome]], h)
  fit0_state <- lm(y_h ~ 0 + state_dummy +
                     shock_high + shock_low + controls_high + controls_low)
  
  nw_vcov_state <- NeweyWest(fit0_state, lag = hac_lag, prewhite = FALSE, adjust = TRUE)
  ct <- coeftest(fit0_state, vcov = nw_vcov_state)
  
  shock_row_high <- ct["shock_high", ]
  shock_row_low  <- ct["shock_low", ]
  
  data.frame(
    h = h,
    coef_high = shock_row_high["Estimate"],
    coef_low = shock_row_low["Estimate"],
    se_high = shock_row_high["Std. Error"],
    se_low = shock_row_low["Std. Error"],
    row.names = NULL
  )
}

#test for h=0 response
estimate_statelp_horizon(
  h = 0, outcome = "y",
  state_dummy = state_lag,
  shock_high = shock_high, shock_low = shock_low,
  controls_high = controls_high, controls_low = controls_low,
  data = rzdat
)


# Generate IRFs -----------------------------------------------------------

irf_y_state <- purrr::map_dfr(0:20, ~ estimate_statelp_horizon(
  h = .x, outcome = "y",
  state_dummy = state_lag,
  shock_high = shock_high, shock_low = shock_low,
  controls_high = controls_high, controls_low = controls_low,
  data = rzdat
))

irf_g_state <- purrr::map_dfr(0:20, ~ estimate_statelp_horizon(
  h = .x, outcome = "g",
  state_dummy = state_lag,
  shock_high = shock_high, shock_low = shock_low,
  controls_high = controls_high, controls_low = controls_low,
  data = rzdat
))

irf_y_state
irf_g_state


# Plot SDLP IRFs ----------------------------------------------------------

plot_irf_state <- function(irf_df, hi_color = "steelblue", lo_color = "firebrick") {
  ggplot(irf_df, aes(x = h)) +
    geom_ribbon(aes(ymin = coef_high - 1.96 * se_high,
                    ymax = coef_high + 1.96 * se_high),
                fill = hi_color, alpha = 0.15) +
    geom_ribbon(aes(ymin = coef_low - 1.96 * se_low,
                    ymax = coef_low + 1.96 * se_low),
                fill = lo_color, alpha = 0.15) +
    geom_line(aes(y = coef_high), color = hi_color, linewidth = 1, linetype = "dashed") +
    geom_line(aes(y = coef_low), color = lo_color, linewidth = 1) +
    geom_point(aes(y = coef_low), color = lo_color, shape = 1, size = 1.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      aspect.ratio = 1
    )
}

p_g_state <- plot_irf_state(irf_g_state)
p_y_state <- plot_irf_state(irf_y_state)
p_g_state / p_y_state

#save
ggsave("output/figures/figure5_state.png", p_g_state / p_y_state, width = 6, height = 8, dpi = 150)




# Blanchard-Perotti Robustness --------------------------------------------

bp_controls <- make_lag_matrix(rzdat, vars = c("y", "g"), p = 4)
dim(bp_controls)
colnames(bp_controls)


estimate_horizon_bp <- function(h, outcome, shock, controls, data, hac_lag = 4) {
  if (outcome == shock & h == 0) {
    return(data.frame(h = h, coef = 1, se = 0, row.names = NULL))
  }
  estimate_horizon(h = h, outcome = outcome, shock = shock,
                   controls = controls, data = data, hac_lag = hac_lag)
}


irf_y_bp <- purrr::map_dfr(0:20, ~ estimate_horizon_bp(h = .x, outcome = "y", shock = "g", controls = bp_controls, data = rzdat))
irf_g_bp <- purrr::map_dfr(0:20, ~ estimate_horizon_bp(h = .x, outcome = "g", shock = "g", controls = bp_controls, data = rzdat))

irf_y_bp
irf_g_bp

