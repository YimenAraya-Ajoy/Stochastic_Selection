# Predicting long-term evolutionary outcomes from fluctuating selection estimates

Code to accompany:

> *Predicting long-term evolutionary outcomes from fluctuating selection estimates: revisiting the optimal clutch size in stochastic environments*

---

## Overview

This repository contains the R code used to fit the statistical models and produce all figures for the manuscript. The analysis links empirical estimates of fluctuating selection — obtained from Bayesian random-regression mixed-effects models — to theoretical parameters governing evolution in stochastic environments, including environmental stochasticity (σ²_e), the covariance penalty *c(z, z̄)*, and the stochastic adaptive topography *m̃(z, z̄)*.

The empirical application revisits a classic dataset on annual variation in clutch size and recruitment in great tits (*Parus major*) at Wytham Woods, Oxford, UK (1980–present), greatly expanding the dataset of Boyce & Perrins (1987).

---

## Repository structure

```
.
├── data/
│   └── data.csv               # Individual-level breeding and recruitment data
├── R/
│   ├── Analysis.R             # Data preparation + Bayesian model fitting (modL, modR) via brms
│   ├── Results.R              # Posterior extraction, derived quantities, and Table 2
│   ├── Figure1.R              # Annual fitness functions + adaptive topographies
│   ├── Figure2.R              # Covariance penalty and environmental stochasticity
│   └── Figure3.R              # Pairwise invasibility plot (PIP)
└── README.md
```

---

## Requirements

### R packages

| Package | Version tested | Purpose |
|---------|---------------|---------|
| `brms`  | ≥ 2.21        | Bayesian multilevel model fitting |
| `posterior` | ≥ 1.6    | Posterior draw manipulation |

Install all dependencies with:

```r
install.packages(c("brms", "posterior"))
```

`brms` requires a working Stan installation. Follow the [RStan getting started guide](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) or use the `cmdstanr` backend:

```r
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()
```

---

## Data

The data were collected at Wytham Woods, Oxford, UK. Each row represents a single breeding attempt, including:

| Column | Description |
|--------|-------------|
| `year` | Breeding year |
| `cs` | Clutch size (number of eggs) |
| `layday` | Date of first egg (day of year) |
| `hatchdate` | Hatch date (day of year) |
| `recruit` | Number of offspring recruiting to the breeding population the following year |

> **Access:** The Wytham Woods great tit data are managed by the University of Oxford.

---

## Workflow

Run scripts in the following order:

### 1. Analysis

```r
source("R/Analysis.R")
```

Loads `data/data.csv`, filters to complete cases on hatch date, globally z-scores clutch size and lay date, computes within-year centred versions of both traits, constructs the composite fitness measure (`r2s = surv + recruit`), and derives relative recruitment (`rel_rec`) and relative fledgling success (`rel_fledge`).

Then fits two complementary Bayesian mixed-effects models:

- **`modL`** — log-linear model for absolute recruitment (negative binomial family). Estimates the mean log-fitness surface and temporal variance–covariance components (V_μ0, V_μ1, C_μ01) for deriving σ²_e(z̄) and the covariance penalty *c(z, z̄)*.
- **`modR`** — arithmetic relative fitness model (Gaussian family, no link function). Analyses recruits divided by the annual mean. Provides the relative fitness selection gradients for comparison.

Both models include fixed linear and quadratic effects of clutch size and lay date, their interaction, and population density as a fixed covariate. Random intercepts and slopes for clutch size and lay date vary by year.

Models are fitted with 3 chains × 4000 iterations each (2000 warm-up). Expect several hours of runtime on a modern laptop.

### 2. Extract posterior quantities

```r
source("R/Results.R")
```

Extracts posterior draws for all variance components, computes σ²_e(z̄) and *c(z, z̄)* across a grid of z̄ values, and derives posterior summaries (mean, 50% and 95% credible intervals) for all plotted quantities, as well as point estimates for the four phenotypic optima reported in the paper (geometric mean, arithmetic mean, arithmetic relative fitness, and stochastic equilibrium z̄*). Also produces Table 2.

### 3. Figures

```r
source("R/Figure1.R")   # Annual fitness functions + four adaptive topographies
source("R/Figure2.R")   # Covariance penalty c(z, z̄) and σ²_e(z̄) vs mean clutch size
source("R/Figure3.R")   # Pairwise invasibility plot
```

Figures are saved as PDFs to the path specified in each script. Update `out_file` at the top of each figure script to match your local directory before running.

---

## Key analytical quantities

The code implements the following expressions from the manuscript. Given posterior draws of V_μ0, V_μ1, and C_μ01 from `modL`:

**Temporal variance in log-fitness of phenotype z** (Eq. 7):
```
Var[r_t(z)] = V_μ0 + z² V_μ1 + 2z C_μ01
```

**Environmental stochasticity** (Eq. 8):
```
σ²_e(z̄) = V_μ0 + z̄² V_μ1 + 2z̄ C_μ01 + ½ σ²_z V_μ1
```

**Covariance penalty** (Eq. 10):
```
c(z, z̄) = V_μ0 + z z̄ V_μ1 + (z + z̄) C_μ01
```

**Stochastic evolutionary equilibrium** (Eq. 14):
```
z̄* = -(β₁ - C_μ01) / (2β₂ - V_μ1)
```

---

## Citation

If you use this code, please cite the manuscript:

> [Authors]. (*in preparation*). Predicting long-term evolutionary outcomes from fluctuating selection estimates: revisiting the optimal clutch size in stochastic environments.

And the `brms` package:

> Bürkner, P.-C. (2017). brms: An R package for Bayesian multilevel models using Stan. *Journal of Statistical Software*, 80, 1–28.

---

## License

Code is released under the [MIT License](LICENSE). Data are subject to separate access conditions — see the Data section above.

---
