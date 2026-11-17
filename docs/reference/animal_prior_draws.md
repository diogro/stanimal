# Prior draws for heritabilities and genetic correlations

Samples the priors that
[`lmm_animal`](https://diogro.github.io/stanimal/reference/lmm_animal.md)
places on the genetic and residual covariance matrices and summarises
them as heritabilities and genetic correlations, which do not depend on
the scale of the data.

## Usage

``` r
animal_prior_draws(
  n_traits,
  n = 4000,
  lkj_prior = 2,
  trait_names = NULL,
  beta_prior = NULL,
  correlation_prior = c("lkj", "size", "shared_size")
)
```

## Arguments

- n_traits:

  Number of traits, a whole number of at least 1.

- n:

  Number of prior draws, a whole number of at least 1.

- lkj_prior:

  LKJ shape used in the fit, a finite number above 0. Shapes below 1 put
  mass on correlations near +/-1, and a shape of 0 or less is not a
  distribution.

- trait_names:

  Optional trait names.

- beta_prior:

  Optional shapes of the beta prior that the `parametrization = "h2"`
  models put on the heritability of each trait. When given, the
  heritabilities come from that beta instead of the scaled half-normal
  ratio, because that is the prior the fit uses. The genetic
  correlations are drawn from the LKJ in either case.

- correlation_prior:

  Prior on the genetic correlations, as in
  [`lmm_animal`](https://diogro.github.io/stanimal/reference/lmm_animal.md):
  `"lkj"` (default), or `"size"` and `"shared_size"`, which draw each
  correlation matrix as \\D(ll^T + R)D\\ with \\R \sim\\
  LKJ(`lkj_prior`), size loadings \\l_k = \mu + \sigma z_k\\,
  \\\mu^2/(1+\mu^2) \sim\\ Beta(1, 1), \\\sigma \sim\\ half-normal(0, 1)
  and \\z_k \sim N(0, 1)\\. The two size priors differ only in how `E`
  shares the loading pattern, so they give the same draws here. Both are
  fitted only with `parametrization = "h2"`, so the heritabilities come
  from `beta_prior`, with `c(1, 1)` when it is `NULL`. They need at
  least two traits.

## Value

A long data frame with columns `type` (`"Heritability"` or
`"Genetic correlation"`), `parameter` and `value`.

## Details

Standard deviations are half-normal(0, 1) and correlation matrices are
LKJ(`lkj_prior`), as in the Stan models. Fits with
`parametrization = "variance"` use that prior unless `beta_prior` is
given, in which case that beta is the prior on each heritability.

The residual correlation structure of the fit does not enter these
draws: heritability divides the genetic variance by the phenotypic
variance, which only uses the diagonal of the residual covariance, and
the genetic correlations come from the correlation matrix of `G`. So
`residual_structure` makes no difference here, and there is no argument
for it.
