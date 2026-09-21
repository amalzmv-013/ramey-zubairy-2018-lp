###############################################
# ramey-zubairy-2018-lp
# Author- Amal Varghese
#  Purpose: Data prep
###############################################


# Prelims -----------------------------------------------------------------

rm(list=ls())
library(readxl); library(dplyr)


# Data loading and screening ----------------------------------------------

#Data Loading
rzdat <- read_excel("data/raw/RZDAT.xlsx", sheet = "rzdat")
glimpse(rzdat)

#Screening for observations post 1889
rzdat <- rzdat %>% filter(quarter >= 1889)
dim (rzdat) 
range(rzdat$quarter)
all(diff(rzdat$quarter) ==0.25)


# Construct real government spending and normalized variables --------

rzdat <- rzdat %>%
  mutate(
    rgov = ngov/pgdp, # real gov spending
    # Gordon–Krenn normalization
    y = rgdp/ rgdp_pott6,
    g = rgov / rgdp_pott6
  )

summary(rzdat$y); summary(rzdat$g)


# Construct news shock ----------------------------------------------------

rzdat <- rzdat %>%
  mutate(
    newsy = news /(dplyr::lag(rgdp_pott6,1)* dplyr::lag(pgdp,1))
  )

summary(rzdat$newsy); head(rzdat$newsy, 3)


# Construct state indicators ----------------------------------------------

rzdat <- rzdat %>%
  mutate(
    wwii= quarter >= 1941.5  & quarter <1946, # WW II
    slack = unemp >= 6.5 # unemployment implied slack regimes
  )

table(rzdat$wwii); table(rzdat$slack, useNA = "always")


# Save --------------------------------------------------------------------

dir.create("data/processed", showWarnings = FALSE)
saveRDS(rzdat, "data/processed/rzdat_prepped.rds")

