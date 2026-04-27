# Short Interest and Aggregate Stock Returns

Replication of Rapach, Ringgenberg, and Zhou (2016), "Short interest and aggregate stock returns," *Journal of Financial Economics*, Vol. 121, pp. 46–65.

## What is replicated

| Output | Script | Description |
|--------|--------|-------------|
| Table 3 | `table_3.R` | In-sample predictive regressions for 14 GW predictors and SII at h = 1, 3, 6, 12 month horizons |
| Table 4 | `table_4.R` | Alternative detrending methods (linear, quadratic, cubic, stochastic) for SII |
| Table 5 | `table_5.R` | Out-of-sample R²OS, Clark–West MSFE-adjusted statistics, and HLN encompassing tests |
| Table 6 | `table_6.R` | Out-of-sample CER gains for a mean-variance investor, full period and subperiods |

All in-sample results cover 1973:01–2014:12. Out-of-sample evaluation period is 1990:01–2014:12.

## How to run

All scripts are R and must be run from the **repo root**.

```r
source("table_3.R")   # Table 3: in-sample predictive regressions + wild bootstrap (~10 min)
source("table_4.R")   # Table 4: alternative detrending
source("table_5.R")   # Table 5: out-of-sample R²OS and encompassing tests
source("table_6.R")   # Table 6: out-of-sample CER gains
```

Each table script sources `load.R`, which reads the data and constructs all variables.

## Methods

### In-sample predictive regressions (Table 3)

Bivariate OLS predictive regressions of the S&P 500 log excess return on each of the 14 Goyal–Welch (2008) predictor variables and SII. Also reports a multiple regression including the first three principal components of the GW predictors and SII. Newey–West HAC *t*-statistics with *h* lags. Wild bootstrapped *p*-values (B = 1,000 fixed-regressor iterations).

### Alternative detrending (Table 4)

SII is constructed by detrending log(EWSI) using four methods — linear, quadratic, cubic, and stochastic (60-month backward-looking moving average) — and standardising to unit variance.

### Out-of-sample forecasts (Tables 5–6)

Expanding-window predictive regression forecasts. SII is recursively detrended to avoid look-ahead bias.

| Metric | Description |
|--------|-------------|
| R²OS | Campbell & Thompson (2008) out-of-sample R²; tested with Clark–West (2007) statistic |
| Encompassing | Harvey–Leybourne–Newbold (1998) test of whether SII-based forecasts encompass GW predictor forecasts |
| CER gain | Annualised certainty-equivalent return gain (basis points) for a mean-variance investor |

### Key parameters

| Parameter | Value | Meaning |
|-----------|-------|---------|
| Sample | 1973:01–2014:12 | 504 monthly observations |
| In-sample end | 1989:12 | R = 204 months |
| OOS evaluation | 1990:01–2014:12 | P = 300 months |
| Forecast horizons | 1, 3, 6, 12 | Monthly, quarterly, semi-annual, annual |
| NW-HAC lags | *h* | Forecast horizon |
| Bootstrap iterations | 1,000 | Fixed-regressor wild bootstrap |
| Risk aversion (γ) | 3 | CER gain calculation |
| Weight bounds | [−0.5, 1.5] | Portfolio weight constraints |
| Volatility window | 120 months | 10-year rolling variance forecast |

## Repository structure

```
data-raw/
  Returns_short_interest_data.xlsx  <- Goyal-Welch + short interest data (1871:01–2014:12)
R/
  nwest-ols.R                       <- nwest_ols(): OLS with Newey-West HAC standard errors
  ols-fit.R                         <- ols_fit(): plain OLS
  constrain-weight.R                <- constrain_weight(): portfolio weight bounds
  compute-cer-sharpe.R              <- compute_cer_sharpe(): CER gains and Sharpe ratios
refs/
  *.pdf                             <- paper (Rapach, Ringgenberg, and Zhou, 2016)
load.R                              <- data preparation; reads XLS, constructs GW predictors, SII, PCs
table_3.R
table_4.R
table_5.R
table_6.R
```

## References

- Rapach, D. E., Ringgenberg, M. C., and Zhou, G. (2016). Short interest and aggregate stock returns. *Journal of Financial Economics*, 121(1), 46–65.
