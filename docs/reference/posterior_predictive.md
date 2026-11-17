# Replicate the response from the posterior

Draws replicate datasets from the fitted model, one per posterior draw,
so that the data can be compared with what the model says they should
look like. This is the individual-level check that
[`fit_loo`](https://diogro.github.io/stanimal/reference/fit_loo.md) is
not: LOO compares model structures by held-out predictive density, while
this shows *which* features of the data the fit reproduces.

## Usage

``` r
posterior_predictive(fit, ndraws = 100, seed = NULL, ...)

# S3 method for class 'stanimal_yrep'
print(x, ...)

# S3 method for class 'stanimal_yrep'
plot(x, ...)
```

## Arguments

- fit:

  A fit from
  [`lmm_animal`](https://diogro.github.io/stanimal/reference/lmm_animal.md)
  or
  [`lm_phen_cov`](https://diogro.github.io/stanimal/reference/lm_phen_cov.md).

- ndraws:

  Number of posterior draws to replicate from; taken evenly across the
  kept draws, so both warmup-free chains are represented.

- seed:

  Optional seed for the simulation, to make the replicates reproducible.
  The sampler's seed does not cover this step.

- ...:

  Unused; for method compatibility.

- x:

  An object returned by `posterior_predictive()`.

## Value

An object of class `stanimal_yrep`: a list with `yrep` (individuals by
traits by draws, the same layout as
[`fit_draws`](https://diogro.github.io/stanimal/reference/fit_draws.md),
with draws last), the observed `y` on the original response scale, the
trait names, and the number of replicates.
[`print()`](https://rdrr.io/r/base/print.html) reports, for the trait
means, standard deviations and pairwise correlations, the observed
value, the 95% interval of the replicates, and `p`: the share of
replicates at or above the observed value, a posterior predictive tail
probability. These probabilities are not generally uniform under a
correct model because the data are used both to fit and to check it.
Values near 0 or 1 can flag a discrepancy; a middling value does not
establish a good fit.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) shows the
replicate distributions.

## Details

For an animal fit the replicate is \\X\beta + u + e\\ with \\u \sim N(0,
A \otimes G)\\ and \\e \sim N(0, I \otimes E)\\; for a phenotypic fit it
is \\X\beta + e\\ with \\e \sim N(0, P)\\. `G`, `E` and `P` are the
reported draws, so they are already on the response scale.

Replicates are simulated for every trait of every individual, including
the entries that were not observed, so `y` is the response as it was
passed in and has `NA` where the data are incomplete. Printed and
plotted summaries compare each replicate on the same observed entries as
the corresponding statistic in `y`. Pairs with fewer than two jointly
observed values or constant observed values are omitted with a warning.

## See also

[`fit_loo`](https://diogro.github.io/stanimal/reference/fit_loo.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pp <- posterior_predictive(fit)
pp                      # observed against replicated summaries
plot(pp)                # and the distributions behind them
} # }
```
