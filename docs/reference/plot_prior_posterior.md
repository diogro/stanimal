# Posterior versus prior for heritabilities and genetic correlations

Posterior versus prior for heritabilities and genetic correlations

## Usage

``` r
plot_prior_posterior(
  fit,
  trait_names = NULL,
  lkj_prior = NULL,
  n_prior = 4000,
  prior = NULL,
  ncol = 4
)

plot_forest(
  fit,
  trait_names = NULL,
  lkj_prior = NULL,
  n_prior = 4000,
  prior = NULL,
  prob = 0.95
)
```

## Arguments

- fit:

  A fit returned by `lmm_animal`, or a `posterior` draws object
  containing `h2` and, for several traits, `corrG`.

- trait_names:

  Optional trait names. Defaults to names recorded on the fit, or
  generic names when these are unavailable.

- lkj_prior:

  Optional LKJ shape for the prior draws. The default `NULL` uses the
  shape recorded by the fit. Supplying a different shape warns that the
  comparison no longer uses the fitted prior. It is ignored, with a
  warning, for factor-G or univariate fits, or when `prior` is given.

- n_prior:

  Passed to
  [`animal_prior_draws`](https://diogro.github.io/stanimal/reference/animal_prior_draws.md)
  as `n`.

- prior:

  Optional prior draws in the format returned by `animal_prior_draws`;
  computed if `NULL` and the fit's prior is known. Its parameter names
  and types must match the posterior. The caller is responsible for
  ensuring the draws represent the fitted prior.

- ncol:

  Number of facet columns for `plot_prior_posterior`.

- prob:

  Probability mass of the highest density interval, between 0 and 1.

## Value

A ggplot object. `plot_prior_posterior` shows overlaid densities in one
panel per parameter, and `plot_forest` shows means and `prob` intervals
of the posterior and the prior.

## Details

Where prior and posterior overlap, the data add little information.

The prior is drawn from the model configuration recorded by
[`lmm_animal`](https://diogro.github.io/stanimal/reference/lmm_animal.md),
so an `"h2"` fit is compared with its beta prior on the heritability and
an unstructured fit with the LKJ shape it was given. A fit with
`correlation_prior = "size"` or `"shared_size"` is compared with the
size-factor prior on its genetic correlations (see
[`animal_prior_draws`](https://diogro.github.io/stanimal/reference/animal_prior_draws.md));
fits made before that argument existed are read as LKJ fits. For
ordinary factor `G`, it simulates the regularized horseshoe loadings
using the recorded factor count, global scale, slab and column decay,
then constructs the genetic correlations. The current model uses the
recorded beta prior on trait heritability; older factor fits recorded as
`"variance"` use independent half-normal genetic and residual scales.
These prior draws use R's random number generator.

The shared-factor model and a fit sampled from a custom `model_file`
cannot be reconstructed here. Bare draws and older fits missing
configuration metadata cannot establish which prior was used either. In
each case pass `prior` with draws from the prior you used, using the
same trait names as the posterior. The same rules apply to
`plot_forest()`.
