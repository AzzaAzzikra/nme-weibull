# nmeweibull

nmeweibull is an R package for modelling positive survival-time data with
the New Modified Exponential-Weibull (NME-Weibull) distribution of
Alshanbari et al. (2023). Its CDF is

$$F(x) = 1 - \frac{\alpha^2 e^{-\lambda x^\beta}}{\left(\alpha + 1 - e^{-\lambda x^\beta}\right)^2}, \qquad x \ge 0,$$

with $\alpha > 0$, shape $\beta > 0$, and rate $\lambda > 0$. The Weibull
distribution $G(x) = 1 - e^{-\lambda x^\beta}$ is the limit as
$\alpha \to \infty$. In terms of `stats::pweibull()`, the scale is
$\sigma = \lambda^{-1/\beta}$.

| Function | Purpose |
|---|---|
| `dnmeweibull()` | Density (PDF) |
| `pnmeweibull()` | Distribution function (CDF), with `lower.tail` and `log.p` |
| `snmeweibull()` | Survival function |
| `hnmeweibull()` | Hazard function |
| `qnmeweibull()` | Quantile function (closed form) |
| `rnmeweibull()` | Random generation (inverse transform) |
| `fitnmeweibull()` | Maximum likelihood estimation for complete data |
| `varnmeweibull()` | Value-at-Risk |
| `tvarnmeweibull()` | Tail Value-at-Risk |
| `rvarnmeweibull()` | Range Value-at-Risk |

## Installation

From the package folder (the one containing `DESCRIPTION`):

```r
install.packages("devtools")
devtools::install_local("path/to/nmeweibull")
library(nmeweibull)
```

## Example

```r
library(nmeweibull)

# Distribution functions
dnmeweibull(1, alpha = 1, beta = 2, lambda = 0.5)
pnmeweibull(1, alpha = 1, beta = 2, lambda = 0.5)
qnmeweibull(0.5, alpha = 1, beta = 2, lambda = 0.5)

# Estimation
set.seed(123)
x <- rnmeweibull(300, alpha = 1, beta = 2, lambda = 0.5)
fit <- fitnmeweibull(x)
fit

# Risk measures from the fitted model
est <- fit$estimate
varnmeweibull(c(0.90, 0.95, 0.99), est["alpha"], est["beta"], est["lambda"])
tvarnmeweibull(c(0.90, 0.95, 0.99), est["alpha"], est["beta"], est["lambda"])
rvarnmeweibull(c(0.90, 0.95), c(0.95, 0.99), est["alpha"], est["beta"],
               est["lambda"])
```

`fitnmeweibull()` reports a convergence code and a `boundary` flag. An
estimate of `alpha` at its upper bound (`1e4` by default) means that the data
cannot be distinguished from a Weibull sample.

## Analysis scripts

The folder `scripts/` holds the analyses reported in the thesis. They are
not part of the installed package. Run them from the package root (open
`nmeweibull.Rproj`); each script writes tables to `results/` and figures to
`figures/`.

| Script | Content |
|---|---|
| `01_distribution_figures.R` | Example values and shapes of the PDF, survival, and hazard functions |
| `02_validation.R` | Validation rules V1-V9 on 48 parameter combinations |
| `03_simulation_mle.R` | Monte Carlo simulation of the MLE (3 scenarios x n = 100, 500, 1000, M = 1000) |
| `04_survival_hazard_gbsg2.R` | Survival and hazard modelling of the GBSG2 recurrence-free survival times, NME-Weibull versus Weibull |
| `05_risk_measures_medical_cost.R` | Actuarial application: model comparison and VaR, TVaR, and RVaR of medical costs |
| `helpers.R` | Shared helper functions used by the scripts |

The data used by the scripts are in `data/`:

- `gbsg2_recurrence_free_time.csv`: recurrence-free survival time in days
  (`time`) and status (`status`: 1 = recurrence or death, 0 = censored) of
  the 686 patients of the German Breast Cancer Study Group 2 trial. Because
  `fitnmeweibull()` handles complete data, the script uses the 299 patients
  with an observed event.
- `medical_cost_charges.csv`: individual medical costs billed by health
  insurance for 1,338 policyholders in US dollars (`charges`), Medical Cost
  Personal Datasets.

## Reference

Alshanbari, H. M., Odhah, O. H., Ahmad, Z., Khan, F., and El-Bagoury,
A. A.-A. H. (2023). A new probability distribution: Model, theory and
analyzing the recovery time data. *Axioms*, 12(5), 477.
<https://doi.org/10.3390/axioms12050477>
