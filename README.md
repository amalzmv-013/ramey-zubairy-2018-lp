# Ramey and Zubairy (2018) Replication

Replication of **Figure 5** and **Figure 6** from:

> Ramey, V. A., & Zubairy, S. (2018). Government spending multipliers in good times and in bad: Evidence from US historical data. *Journal of Political Economy*, 126(2), 850–901.

## What this replicates

Impulse response functions from **linear and state-dependent local projections**, estimating the response of government spending and GDP to a military news shock. States are defined by the unemployment rate threshold (6.5%), following RZ (2018).

## Data

`data/raw/RZDAT.xlsx` — Ramey's publicly available dataset (revised Nov 2016), obtained from the UCSD website.

## How to run

Scripts are in `R/` and should be run in numbered order:

| Script | Contents |
|------------------------------------|------------------------------------|
| `01_data_prep.R` | Load Excel, construct real series, scale by trend GDP |
| `02_linear_lp.R` | Linear LP estimation, Newey-West SEs |
| `03_state_dep_lp.R` | State-dependent LP, unemployment threshold, WWII exclusion |
| `04_figures.R` | Reproduce Figure 5 panels |

## Packages required

``` r
install.packages(c("readxl", "sandwich", "ggplot2", "patchwork", "dplyr"))
```

## Results

Replicated figures are in `output/figures/`.
