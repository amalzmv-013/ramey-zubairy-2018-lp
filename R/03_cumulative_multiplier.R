###############################################
# ramey-zubairy-2018-lp
# Author- Amal Varghese
#  Purpose: Cumulative Multiplier
###############################################


# Prelims -----------------------------------------------------------------

rm(list = ls())
library(dplyr); library(sandwich); library(lmtest)

rzdat<- readRDS("data/processed/rzdat_prepped.rds")
dim(rzdat)


# Forward Cumulative Sums -------------------------------------------------

cum_forward <- function(x,h){
  rowSums(sapply(0:h, function(j) dplyr::lead(x,j)))
}


# Controls ----------------------------------------------------------------

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



# cumulative variables at the 2-year horizon ------------------------------

h<-7
cumy <- cum_forward(rzdat$y, h)
cumg <- cum_forward(rzdat$g, h)

keep <- complete.cases(cumy, cumg, rzdat$newsy, lag_controls)
sum(keep)
which(!keep)




# 2SLS vs ivreg -----------------------------------------------------------

#2SLS

cumy_est     <- cumy[keep]
cumg_est     <- cumg[keep]
shock_est    <- rzdat$newsy[keep]
controls_est <- lag_controls[keep, ]

first_stage  <- lm(cumg_est ~ shock_est + controls_est)
reduced_form <- lm(cumy_est ~ shock_est + controls_est)

cumg_hat     <- fitted(first_stage)
second_stage <- lm(cumy_est ~ cumg_hat + controls_est)

coef(reduced_form)["shock_est"] / coef(first_stage)["shock_est"]
coef(second_stage)["cumg_hat"]


#ivreg
iv_fit <- AER::ivreg(cumy_est ~ cumg_est + controls_est | shock_est + controls_est)

coef(iv_fit)["cumg_est"]
summary(second_stage)$coefficients["cumg_hat", ]
summary(iv_fit)$coefficients["cumg_est", ]



# HAC std error -----------------------------------------------------------

bw_auto <- sandwich::bwNeweyWest(iv_fit, prewhite = FALSE)
bw_auto

nw_iv <- sandwich::NeweyWest(iv_fit, lag = NULL, prewhite = FALSE, adjust = FALSE)
ct_iv <- lmtest::coeftest(iv_fit, vcov = nw_iv)
ct_iv["cumg_est", ]


# F test ------------------------------------------------------------------

nw_fs <- sandwich::NeweyWest(first_stage, lag = NULL, prewhite = FALSE, adjust = FALSE)
ct_fs <- lmtest::coeftest(first_stage, vcov = nw_fs)
F_robust <- ct_fs["shock_est", "t value"]^2
F_robust


lags_to_try <- c(0, 4, 8, 12, 16)

se_by_lag <- sapply(lags_to_try, function(L) {
  v <- sandwich::NeweyWest(iv_fit, lag = L, prewhite = FALSE, adjust = FALSE)
  sqrt(diag(v))["cumg_est"]
})
names(se_by_lag) <- paste0("L", lags_to_try)
se_by_lag




# wrap one horizon into a function ----------------------------------------

estimate_cum_multiplier <- function(h, data, shock, controls, hac_lag = NULL) {
  cumy <- cum_forward(data$y, h)
  cumg <- cum_forward(data$g, h)
  shock_var <- data[[shock]]
  keep <- complete.cases(cumy, cumg, shock_var, controls)
  
  cumy_est     <- cumy[keep]
  cumg_est     <- cumg[keep]
  shock_est    <- shock_var[keep]
  controls_est <- controls[keep, ]
  
  iv_fit <- AER::ivreg(cumy_est ~ cumg_est + controls_est | shock_est + controls_est)
  nw_iv  <- sandwich::NeweyWest(iv_fit, lag = hac_lag, prewhite = FALSE, adjust = FALSE)
  ct     <- lmtest::coeftest(iv_fit, vcov = nw_iv)
  
  first_stage <- lm(cumg_est ~ shock_est + controls_est)
  nw_fs       <- sandwich::NeweyWest(first_stage, lag = hac_lag, prewhite = FALSE, adjust = FALSE)
  ct_fs       <- lmtest::coeftest(first_stage, vcov = nw_fs)
  
  data.frame(
    h        = h,
    mult     = ct["cumg_est", "Estimate"],
    se       = ct["cumg_est", "Std. Error"],
    F_robust = ct_fs["shock_est", "t value"]^2,
    n        = sum(keep),
    row.names = NULL
  )
}

#testing
estimate_cum_multiplier(h = 7,  data = rzdat, shock = "newsy", controls = lag_controls)
estimate_cum_multiplier(h = 15, data = rzdat, shock = "newsy", controls = lag_controls)

#looping over horizons 

irf_mult_linear <- purrr::map_dfr(0:20, ~ estimate_cum_multiplier(
  h = .x, data = rzdat, shock = "newsy", controls = lag_controls
))
irf_mult_linear


# Plot --------------------------------------------------------------------

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


plot_F <- function(df, threshold = 23.1, cap = 30) {
  ggplot(df, aes(x = h, y = pmin(F_robust - threshold, cap))) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          aspect.ratio = 1)
}

p_mult_linear <- plot_irf(dplyr::rename(irf_mult_linear, coef = mult))
p_F_linear    <- plot_F(irf_mult_linear)

p_mult_linear
p_F_linear




# validate the one-step IV against the three-step ratio -------------------


H <- 7
keep_H <- complete.cases(cum_forward(rzdat$y, H), cum_forward(rzdat$g, H),
                         rzdat$newsy, lag_controls)
shock_H    <- rzdat$newsy[keep_H]
controls_H <- lag_controls[keep_H, ]

beta_on_shock <- function(x, j, keep, shock, controls) {
  out_j <- dplyr::lead(x, j)[keep]
  coef(lm(out_j ~ shock + controls))["shock"]
}

sum_beta_y <- sum(sapply(0:H, function(j) beta_on_shock(rzdat$y, j, keep_H, shock_H, controls_H)))
sum_beta_g <- sum(sapply(0:H, function(j) beta_on_shock(rzdat$g, j, keep_H, shock_H, controls_H)))

sum_beta_y / sum_beta_g


# Save --------------------------------------------------------------------

ggsave("output/figures/figure6_linear.png", p_mult_linear, width = 4, height = 4, dpi = 150)
ggsave("output/figures/figure4_F_linear.png", p_F_linear, width = 4, height = 4, dpi = 150)



# State wise results ------------------------------------------------------

state_lag <- dplyr::lag(rzdat$slack, 1)
state_num <- as.numeric(state_lag)

shock_high    <- rzdat$newsy * state_num
shock_low     <- rzdat$newsy * (1 - state_num)
controls_high <- lag_controls * state_num
controls_low  <- lag_controls * (1 - state_num)

colnames(controls_high) <- paste0(colnames(lag_controls), "_high")
colnames(controls_low)  <- paste0(colnames(lag_controls), "_low")

exog_state <- cbind(state_num, controls_high, controls_low)

table(state_num, useNA = "always")
dim(exog_state)


# Multipliers

h <- 7
cumy <- cum_forward(rzdat$y, h)
cumg <- cum_forward(rzdat$g, h)
cumg_high <- cumg * state_num
cumg_low  <- cumg * (1 - state_num)

keep <- complete.cases(cumy, cumg, shock_high, shock_low, exog_state)
sum(keep)

cumy_k       <- cumy[keep]
cumg_high_k  <- cumg_high[keep]
cumg_low_k   <- cumg_low[keep]
shock_high_k <- shock_high[keep]
shock_low_k  <- shock_low[keep]
exog_k       <- exog_state[keep, ]

iv_high <- AER::ivreg(cumy_k ~ cumg_high_k + exog_k | shock_high_k + exog_k)
iv_low  <- AER::ivreg(cumy_k ~ cumg_low_k  + exog_k | shock_low_k  + exog_k)

hac_ct <- function(fit) {
  lmtest::coeftest(fit, vcov = sandwich::NeweyWest(fit, lag = NULL, prewhite = FALSE, adjust = FALSE))
}

hac_ct(iv_high)["cumg_high_k", ]
hac_ct(iv_low)["cumg_low_k", ]

#the joint regression and the test of equal multipliers
iv_joint <- AER::ivreg(cumy_k ~ cumg_high_k + cumg_low_k + exog_k |
                         shock_high_k + shock_low_k + exog_k)
V_joint  <- sandwich::NeweyWest(iv_joint, lag = NULL, prewhite = FALSE, adjust = FALSE)

coef(iv_joint)[c("cumg_high_k", "cumg_low_k")]

diff_hat <- coef(iv_joint)["cumg_high_k"] - coef(iv_joint)["cumg_low_k"]
var_diff <- V_joint["cumg_high_k", "cumg_high_k"] +
  V_joint["cumg_low_k",  "cumg_low_k"]  -
  2 * V_joint["cumg_high_k", "cumg_low_k"]
z_diff   <- diff_hat / sqrt(var_diff)
p_diff   <- 2 * pnorm(-abs(z_diff))

c(diff = unname(diff_hat), se = sqrt(var_diff), p = unname(p_diff))


# function definition

estimate_cum_multiplier_state <- function(h, data, state, shock_high, shock_low,
                                          exog, hac_lag = NULL) {
  cumy      <- cum_forward(data$y, h)
  cumg      <- cum_forward(data$g, h)
  cumg_high <- cumg * state
  cumg_low  <- cumg * (1 - state)
  keep <- complete.cases(cumy, cumg, shock_high, shock_low, exog)
  
  cumy_k       <- cumy[keep]
  cumg_high_k  <- cumg_high[keep]
  cumg_low_k   <- cumg_low[keep]
  shock_high_k <- shock_high[keep]
  shock_low_k  <- shock_low[keep]
  exog_k       <- exog[keep, ]
  
  hac <- function(fit) {
    sandwich::NeweyWest(fit, lag = hac_lag, prewhite = FALSE, adjust = FALSE)
  }
  
  iv_high  <- AER::ivreg(cumy_k ~ cumg_high_k + exog_k | shock_high_k + exog_k)
  iv_low   <- AER::ivreg(cumy_k ~ cumg_low_k  + exog_k | shock_low_k  + exog_k)
  iv_joint <- AER::ivreg(cumy_k ~ cumg_high_k + cumg_low_k + exog_k |
                           shock_high_k + shock_low_k + exog_k)
  
  ct_high <- lmtest::coeftest(iv_high, vcov = hac(iv_high))
  ct_low  <- lmtest::coeftest(iv_low,  vcov = hac(iv_low))
  
  V        <- hac(iv_joint)
  b        <- coef(iv_joint)
  diff_hat <- b["cumg_high_k"] - b["cumg_low_k"]
  var_diff <- V["cumg_high_k", "cumg_high_k"] + V["cumg_low_k", "cumg_low_k"] -
    2 * V["cumg_high_k", "cumg_low_k"]
  
  fs_high    <- lm(cumg_high_k ~ shock_high_k + exog_k)
  fs_low     <- lm(cumg_low_k  ~ shock_low_k  + exog_k)
  V_fs_high <- hac(fs_high)
  V_fs_low  <- hac(fs_low)
  F_high <- unname(coef(fs_high)["shock_high_k"]^2 / V_fs_high["shock_high_k", "shock_high_k"])
  F_low  <- unname(coef(fs_low)["shock_low_k"]^2  / V_fs_low["shock_low_k",  "shock_low_k"])
  
  data.frame(
    h         = h,
    mult_high = ct_high["cumg_high_k", "Estimate"],
    se_high   = ct_high["cumg_high_k", "Std. Error"],
    mult_low  = ct_low["cumg_low_k", "Estimate"],
    se_low    = ct_low["cumg_low_k", "Std. Error"],
    F_high    = F_high,
    F_low     = F_low,
    p_diff    = 2 * pnorm(-abs(unname(diff_hat) / sqrt(var_diff))),
    n         = sum(keep),
    row.names = NULL
  )
}


#testing
estimate_cum_multiplier_state(h = 7,  data = rzdat, state = state_num,
                              shock_high = shock_high, shock_low = shock_low, exog = exog_state)
estimate_cum_multiplier_state(h = 15, data = rzdat, state = state_num,
                              shock_high = shock_high, shock_low = shock_low, exog = exog_state)

#Loop

irf_mult_state <- purrr::map_dfr(0:20, ~ estimate_cum_multiplier_state(
  h = .x, data = rzdat, state = state_num,
  shock_high = shock_high, shock_low = shock_low, exog = exog_state
))
irf_mult_state


# Plot (Fig 6)


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

p_mult_state <- plot_irf_state(
  dplyr::rename(irf_mult_state,
                coef_high = mult_high, coef_low = mult_low)
) + coord_cartesian(ylim = c(-5, 2))

p_mult_state

# F stat plot (Fig 4)


p_F_high <- plot_F(dplyr::rename(irf_mult_state, F_robust = F_high))
p_F_low  <- plot_F(dplyr::rename(irf_mult_state, F_robust = F_low))

zoom_F <- coord_cartesian(ylim = c(-25, 40))

fig4_news <- (p_F_linear + zoom_F + labs(x = "h") + ggtitle("Linear")) |
  (p_F_high   + zoom_F + labs(x = "h") + ggtitle("High unemployment")) |
  (p_F_low    + zoom_F + labs(x = "h") + ggtitle("Low unemployment"))

fig4_news

#save
ggsave("output/figures/figure4_F_news.png", fig4_news, width = 9, height = 3.5, dpi = 150)




# Figure 6  ---------------------------------------------------------------

shift_q <- function(df) dplyr::mutate(df, h = h + 1)

p_top <- plot_irf(dplyr::rename(shift_q(irf_mult_linear), coef = mult)) +
  ggtitle("Linear: cumulative spending multiplier")

p_bot <- plot_irf_state(dplyr::rename(shift_q(irf_mult_state),
                                      coef_high = mult_high, coef_low = mult_low)) +
  coord_cartesian(ylim = c(-5, 2)) +
  ggtitle("State dependent: cumulative spending multiplier") +
  labs(x = "quarter")

fig6 <- p_top / p_bot
fig6

ggsave("output/figures/figure6_multipliers.png", fig6, width = 5, height = 8, dpi = 150)







# NBER State  -------------------------------------------------------------

class(rzdat$recession)
table(rzdat$recession, useNA = "always")
table(dplyr::lag(rzdat$recession, 1), useNA = "always")




make_state_inputs <- function(state, shock, controls) {
  controls_high <- controls * state
  controls_low  <- controls * (1 - state)
  colnames(controls_high) <- paste0(colnames(controls), "_high")
  colnames(controls_low)  <- paste0(colnames(controls), "_low")
  list(
    state      = state,
    shock_high = shock * state,
    shock_low  = shock * (1 - state),
    exog       = cbind(state_num = state, controls_high, controls_low)
  )
}

inputs_unemp <- make_state_inputs(state_num, rzdat$newsy, lag_controls)
all.equal(inputs_unemp$exog, exog_state)

state_nber  <- as.numeric(dplyr::lag(rzdat$recession, 1))
inputs_nber <- make_state_inputs(state_nber, rzdat$newsy, lag_controls)

est_nber <- function(h) {
  estimate_cum_multiplier_state(
    h = h, data = rzdat, state = inputs_nber$state,
    shock_high = inputs_nber$shock_high, shock_low = inputs_nber$shock_low,
    exog = inputs_nber$exog)
}
est_nber(7)
est_nber(15)




# replication-versus-paper table ------------------------------------------

est_state <- function(h, inputs) {
  estimate_cum_multiplier_state(
    h = h, data = rzdat, state = inputs$state,
    shock_high = inputs$shock_high, shock_low = inputs$shock_low,
    exog = inputs$exog)
}

lin7  <- estimate_cum_multiplier(7,  rzdat, "newsy", lag_controls)
lin15 <- estimate_cum_multiplier(15, rzdat, "newsy", lag_controls)
un7   <- est_state(7,  inputs_unemp);  un15 <- est_state(15, inputs_unemp)
nb7   <- est_state(7,  inputs_nber);   nb15 <- est_state(15, inputs_nber)

tab_row <- function(spec, horizon, ours, ours_se, paper, paper_se = NA) {
  data.frame(spec, horizon, ours = round(ours, 3), ours_se = round(ours_se, 3),
             paper, paper_se)
}

comparison <- rbind(
  tab_row("Linear",             "2-year", lin7$mult,      lin7$se,      .66, .067),
  tab_row("Linear",             "4-year", lin15$mult,     lin15$se,     .71, .044),
  tab_row("Unemployment: high", "2-year", un7$mult_high,  un7$se_high,  .60, .095),
  tab_row("Unemployment: low",  "2-year", un7$mult_low,   un7$se_low,   .59, .091),
  tab_row("Unemployment: high", "4-year", un15$mult_high, un15$se_high, .68, .052),
  tab_row("Unemployment: low",  "4-year", un15$mult_low,  un15$se_low,  .67, .121),
  tab_row("NBER: recession",    "2-year", nb7$mult_high,  nb7$se_high,  .63),
  tab_row("NBER: expansion",    "2-year", nb7$mult_low,   nb7$se_low,   .55),
  tab_row("NBER: recession",    "4-year", nb15$mult_high, nb15$se_high, .67),
  tab_row("NBER: expansion",    "4-year", nb15$mult_low,  nb15$se_low,  .64)
)
comparison

dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(comparison, "output/tables/multiplier_comparison.csv", row.names = FALSE)
