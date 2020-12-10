
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bicopula

<!-- badges: start -->

[![R build
status](https://github.com/ellessenne/bicopula/workflows/R-CMD-check/badge.svg)](https://github.com/ellessenne/bicopula/actions)
[![Codecov test
coverage](https://codecov.io/gh/ellessenne/bicopula/branch/master/graph/badge.svg)](https://codecov.io/gh/ellessenne/bicopula?branch=master)<!-- badges: end -->

## Installation

You can install the development version of {bicopula} from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ellessenne/bicopula")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(bicopula)

set.seed(42)
N <- 200
df <- data.frame(
  age = runif(N, 20, 40),
  sex = rbinom(N, size = 1, prob = 0.5),
  biomarker = runif(N, 0, 120)
)

simdata <- simulate_bisurv(
  dist1 = Exponential$new(lambda = 0.1),
  dist2 = Weibull$new(lambda = 0.1, gamma = 1.5),
  formula1 = ~ age + sex,
  formula2 = ~ age + biomarker,
  data = df,
  beta1 = c(age = 0.01, sex = 0.5),
  beta2 = c(age = 0.01, biomarker = 0.005),
  copula = copula::frankCopula(param = 5)
)

with(simdata, plot(eventtime1, eventtime2))
```

<img src="man/figures/README-example-1.png" width="90%" />
