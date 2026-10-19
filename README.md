# glt

**Bayesian Factor Analysis via the Generalized Lower Triangular Process Priors**

## Overview

The package `glt` implements sparse Bayesian factor analysis under the **generalized lower triangular process prior** on the loading matrix by [Canale and Frühwirth-Schnatter, (2026)](https://arxiv.org/abs/2610.09762).

The number of active factors, their sparsity pattern, and the loadings themselves are learned jointly via a Markov chain Monte Carlo sampler implemented in C++ with `Rcpp` and `RcppArmadillo`.

Install the R package with 

```R
install.packages("devtools")
devtools::install_github("tonycanale/glt")
```

## Author

Antonio Canale


## Details

Currently under internal checks at CRAN. 