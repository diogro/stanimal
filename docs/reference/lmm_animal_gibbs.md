# Fit an animal model with a Gibbs sampler

Fits a univariate or multivariate linear animal model with fixed effects
and relatedness matrix using a single-chain Gibbs sampler. This is
intended for quick simulations and checks; its priors differ from those
in `lmm_animal`.

## Usage

``` r
lmm_animal_gibbs(Y, X = NULL, A, iter = 1000, warmup = floor(iter/2), thin = 1)
```

## Arguments

- Y:

  Numeric response vector or matrix, with individuals in rows. Every
  trait must vary, and is checked as in
  [`lmm_animal`](https://diogro.github.io/stanimal/reference/lmm_animal.md).

- X:

  Optional numeric full-column-rank fixed-effects design matrix, with
  the same rows as `Y`. The default `NULL` fits an intercept-only model:
  one column of ones, named `"(Intercept)"`. Unlike the Stan models, it
  needs no column representing a constant: the sampler does not shift
  the response mean through the design.

- A:

  Symmetric positive-semidefinite relationship matrix, with rows and
  columns in the same individual order as the rows of `Y`. When both `A`
  and `Y` have row names they are used to put `A` in the order of `Y`;
  names that disagree are an error.

- iter:

  Total number of Gibbs iterations, a whole number.

- warmup:

  Number of initial iterations to discard, a whole number smaller than
  `iter`.

- thin:

  Keep every `thin`th iteration after burn-in, a whole number;
  `warmup + thin` must not exceed `iter`.

## Value

A list of posterior draws: `G` and `R` are covariance arrays, `beta`
contains fixed-effect matrix draws, `h2` contains trait-wise variance
ratios, and `iterations` gives the retained iteration numbers.
Quantities are indexed by trait first and draw last, from named
dimensions, the layout that
[`back_transform`](https://diogro.github.io/stanimal/reference/back_transform.md)
returns as well.
[`fit_draws`](https://diogro.github.io/stanimal/reference/fit_draws.md)
accepts this list and returns the same layout with `R` named `E`, plus
the quantities derived from it.

## Details

Fixed effects have an improper flat prior. Genetic and residual
covariance matrices have inverse-Wishart priors with degrees of freedom
`K + 1` and scale matrix `(K + 1) I / 2`, where `K` is the number of
traits. These differ from the priors used by `lmm_animal`; the resulting
posterior draws are not directly interchangeable.
