# Posterior means or quantiles of a fit, as matrices

Summarises the draws of an animal model without the draw dimension, so
each covariance matrix is returned as a trait by trait matrix and each
trait vector as a named vector. These are the shapes used to compare
estimates with a simulated truth, compute heritabilities, or look at a
correlation matrix.

## Usage

``` r
draws_mean(x, parameter = NULL, traits = NULL)

draws_quantile(x, parameter = NULL, probs = c(0.025, 0.975), traits = NULL)
```

## Arguments

- x:

  A fit, or an object from
  [`fit_draws`](https://diogro.github.io/stanimal/reference/fit_draws.md).

- parameter:

  Name of the quantity to summarise: `"G"`, `"E"`, `"P"`, `"corrG"`,
  `"corrE"`, `"corrP"`, `"h2"`, `"sigma_G"`, `"sigma_E"`, `"sigma_P"`,
  `"loadings"`, `"uniqueness"`, `"communality"`, `"genetic_uniqueness"`,
  `"residual_uniqueness"`, `"s_G"`, `"factor_heritability"`,
  `"factor_magnitude"` or `"beta"`. The default summarises every
  quantity the fit holds.

- traits:

  Optional trait names, one per trait. The default uses the column names
  of the response when the fit kept them, and `trait1`, `trait2`, ...
  otherwise.

- probs:

  For `draws_quantile`, one or more probabilities between 0 and 1; the
  default is the 95% equal-tailed interval.

## Value

For one `parameter`, a named matrix or vector: `k` by `k` for a
covariance or correlation matrix, `k` by the number of terms for `beta`,
and a vector named after the traits for the others. For `draws_quantile`
the last dimension is the quantile, named after the probability. When
`parameter` is `NULL`, a named list with one summary per quantity.

## Details

The summaries are computed from the draws, so the posterior mean of `G`
is the average of the draws and not the covariance of any average.
`draws_quantile` returns the same shapes as `draws_mean`, with a final
dimension for the probabilities.

## See also

[`fit_draws`](https://diogro.github.io/stanimal/reference/fit_draws.md)

## Examples

``` r
k <- 2
vars <- sprintf("G[%d,%d]", rep(1:k, k), rep(1:k, each = k))
m <- cbind(matrix(stats::rnorm(40), 10, 4, dimnames = list(NULL, vars)),
           matrix(stats::rnorm(40), 10, 4,
                  dimnames = list(NULL, sub("G", "E", vars))))
draws_mean(m, "G", traits = c("t1", "t2"))
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
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#>             t1          t2
#> t1 -0.45074303 -0.09215854
#> t2 -0.05647442 -0.13190713
round(draws_mean(m, "G", traits = c("t1", "t2")), 2)
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
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#>       t1    t2
#> t1 -0.45 -0.09
#> t2 -0.06 -0.13
draws_quantile(m, "h2", traits = c("t1", "t2"))
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
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#>         2.5%    97.5%
#> t1 -1.562079 1.211529
#> t2 -6.552883 2.804419
```
