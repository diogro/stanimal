# Posterior draws of a fit, shaped by trait

Extracts the draws of an animal model and reshapes every parameter so
that its first dimensions are the traits. Covariance matrices come out
as trait by trait by draw arrays with named dimensions, ready to
average, take correlations of, or compare with a simulated truth,
instead of one Stan variable per element such as `G[1,1]`.

## Usage

``` r
fit_draws(fit, traits = NULL, terms = NULL)
```

## Arguments

- fit:

  A fit from
  [`lmm_animal`](https://diogro.github.io/stanimal/reference/lmm_animal.md),
  [`lm_phen_cov`](https://diogro.github.io/stanimal/reference/lm_phen_cov.md)
  or
  [`lmm_animal_gibbs`](https://diogro.github.io/stanimal/reference/lmm_animal_gibbs.md),
  a `posterior` draws object, or a named matrix of draws such as
  `fit$draws(format = "draws_matrix")`. Objects from
  [`back_transform`](https://diogro.github.io/stanimal/reference/back_transform.md)
  also work.

- traits:

  Optional trait names, one per trait. The default uses the column names
  of the response when the fit kept them, and `trait1`, `trait2`, ...
  otherwise.

- terms:

  Optional names for the columns of `X` in the `beta` draws, in their
  order in the fit. The default uses the column names of `X` when the
  fit kept them, and the names the fit reports otherwise.

## Value

A list of class `fit_draws`, with one element per quantity the fit
holds: `G`, `E`, `P`, `corrG`, `corrE`, `corrP` as trait by trait by
draw arrays, `h2`, `sigma_G`, `sigma_E`, `sigma_P`, `uniqueness`,
`communality`, `genetic_uniqueness`, `residual_uniqueness`, `s_G` and,
for `correlation_prior = "size"` or `"shared_size"` fits,
`size_loadings_G` and `size_loadings_E` as trait by draw matrices,
`loadings` as a trait by factor by draw array for factor fits,
`factor_heritability` and `factor_magnitude` as factor by draw matrices,
and `beta` as a trait by term by draw array. Every dimension is named,
and draws are last, the same layout as
[`lmm_animal_gibbs`](https://diogro.github.io/stanimal/reference/lmm_animal_gibbs.md)
and
[`back_transform`](https://diogro.github.io/stanimal/reference/back_transform.md).
Quantities the fit does not report are derived when they follow from the
others: `P = G + E`, `corrG`, `corrE` and `corrP` from the covariances,
`sigma_G`, `sigma_E` and `sigma_P` from their diagonals, and `h2` from
`G` and `P`.

## Details

The Gibbs sampler calls its residual covariance `R`; it is returned as
`E` so that the two engines give the same names. The priors differ
between the engines, so the draws themselves are not interchangeable.

`loadings` comes back \*\*sign-aligned\*\*: a column of \\\Lambda\\ is
identified only up to its sign, so a chain may visit both orientations,
and averaging the raw draws then cancels the factor out of the posterior
mean and out of any heatmap built from it. Each draw is flipped so that
the trait with the largest mean absolute loading on that factor is
positive, which is a property of the whole posterior rather than of one
draw. Nothing else is changed: `communality`, `uniqueness`,
`factor_magnitude` and the covariances are invariant to the sign anyway.
A permutation of the columns or a genuine rotation between draws is
\*not\* repaired, which is why those are read on `G` and `communality`
instead. For the draws exactly as the sampler wrote them, use
`fit$draws("loadings")`.

## See also

[`draws_mean`](https://diogro.github.io/stanimal/reference/draws_mean.md),
[`draws_quantile`](https://diogro.github.io/stanimal/reference/draws_mean.md),
[`animal_parameter_draws`](https://diogro.github.io/stanimal/reference/animal_parameter_draws.md)

## Examples

``` r
k <- 2
vars <- sprintf("G[%d,%d]", rep(1:k, k), rep(1:k, each = k))
m <- cbind(matrix(stats::rnorm(40), 10, 4, dimnames = list(NULL, vars)),
           matrix(stats::rnorm(40), 10, 4,
                  dimnames = list(NULL, sub("G", "E", vars))))
d <- fit_draws(m, traits = c("t1", "t2"))
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
names(d)
#>  [1] "G"       "E"       "P"       "corrG"   "corrE"   "corrP"   "h2"     
#>  [8] "sigma_G" "sigma_E" "sigma_P"
dim(d$G)                       # traits x traits x draws
#> [1]  2  2 10
round(draws_mean(d, "G"), 2)   # posterior mean G, as a matrix
#>      t1   t2
#> t1 0.07 0.27
#> t2 0.24 0.12
round(draws_mean(d, "corrG"), 2)
#>     t1  t2
#> t1   1 NaN
#> t2 NaN   1
```
