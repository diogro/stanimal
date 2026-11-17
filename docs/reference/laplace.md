# Laplace approximations of the covariance fits

Fits the same model as
[`lmm_animal`](https://diogro.github.io/stanimal/reference/lmm_animal.md)
or
[`lm_phen_cov`](https://diogro.github.io/stanimal/reference/lm_phen_cov.md)
and returns a Gaussian approximation at the posterior mode instead of
posterior draws. It can be much faster than sampling, but has no MCMC
diagnostics: there is no `Rhat`, divergence count or effective sample
size, because there are no chains.

## Usage

``` r
lmm_animal_laplace(
  Y,
  X = NULL,
  A,
  parametrization = c("variance", "h2"),
  lkj_prior = 2,
  beta_prior = c(1, 1),
  blend_A = FALSE,
  model_file = NULL,
  residual_structure = c("unstructured", "diagonal", "shared_factor"),
  genetic_structure = c("unstructured", "factor"),
  n_factors = NULL,
  hs_global_scale = 1,
  hs_slab_scale = 1,
  column_decay = 0.8,
  draws = 2000,
  seed = NULL,
  correlation_prior = c("lkj", "size", "shared_size"),
  ...
)

lm_phen_cov_laplace(
  Y,
  X = NULL,
  lkj_prior = 2,
  covariance_structure = c("unstructured", "factor"),
  n_factors = NULL,
  hs_global_scale = 1,
  hs_slab_scale = 1,
  column_decay = 0.8,
  model_file = NULL,
  draws = 2000,
  seed = NULL,
  ...
)
```

## Arguments

- Y, X, A:

  Relationship matrix and data, as in
  [`lmm_animal`](https://diogro.github.io/stanimal/reference/lmm_animal.md).
  `lmm_animal_laplace()` needs `A`; `lm_phen_cov_laplace()` does not
  take one.

- parametrization:

  Which parameters carry the prior. `"variance"` (default) estimates the
  genetic and residual standard deviations with half-normal priors and
  derives h2. `"h2"` estimates a phenotypic standard deviation and a
  heritability per trait, and puts the beta prior (`beta_prior`) on h2.
  It is available for the univariate models and for the unstructured
  multivariate model, and the two differ only in the prior they imply on
  h2: because the response is standardized trait by trait, \\s_G^2 +
  s_E^2\\ is the phenotypic variance, so the half-normal priors of
  `"variance"` are priors on a circle and they induce \\p(h^2) \propto
  1/\sqrt{h^2(1-h^2)}\\, a density whose mass sits at the two ends. The
  two parameterizations agree when the data are informative and separate
  when they are not, so prefer `"h2"` when there are few individuals
  relative to the number of traits, when h2 is the object of inference
  rather than `G`, or when a `"variance"` fit puts its heritabilities
  against 0 or 1. No sample-size threshold has been validated, so treat
  that as a reason to check rather than a rule. The ordinary factor
  model has no switch: its scales are locked to the phenotypic scale and
  it always puts the beta prior on the trait heritabilities. Supplying
  this argument for the factor or shared-factor model warns that it is
  ignored, and a univariate `"variance"` fit with low E-BFMI already
  warns and suggests `"h2"`.

- lkj_prior:

  Positive shape parameter of the LKJ prior. `1` is uniform over
  correlation matrices, not over each pairwise correlation: their
  marginal distribution concentrates increasingly around zero as
  dimension increases. Larger shapes add shrinkage toward zero. For
  multivariate models it applies to unstructured `G`, including with
  diagonal `E`, and to unstructured `E`, including with factor `G`. It
  is unused for factor `G` with diagonal `E`, shared-factor and
  univariate models; explicitly supplying it there warns.

- beta_prior:

  Parameters of the beta prior on h2. Used by the univariate model with
  `parametrization = "h2"`, by the unstructured multivariate model with
  `parametrization = "h2"`, by the ordinary factor model for the trait
  heritabilities, and by the shared-factor model for the factor
  heritabilities. The default `c(1, 1)` is flat on h2 and on each factor
  heritability. It is validated for every model; explicitly supplying it
  where it is not used warns that it is ignored.

- blend_A:

  If TRUE, replace A with 0.99 \* A + 0.01 \* I before fitting, making
  singular A positive definite.

- model_file:

  Optional path to a custom Stan file replacing the default model.

- residual_structure:

  Residual covariance for multivariate fits: "unstructured" (default,
  full covariance between traits), "diagonal", or "shared_factor". With
  a diagonal residual, environmental covariances between traits are
  absorbed into the genetic covariance and inflate heritability, so
  "unstructured" is recommended when estimating G without a factor
  model. "shared_factor" requires `genetic_structure = "factor"` and
  fits the shared-factor model \\G = \Lambda
  \mathrm{diag}(h^2\_\text{factor}) \Lambda^T + \Psi_G\\, \\E = \Lambda
  \mathrm{diag}(1 - h^2\_\text{factor}) \Lambda^T + \Psi_E\\, so the
  same latent factors explain genetic and residual covariance and each
  factor has a heritability.

- genetic_structure:

  Genetic covariance for multivariate fits: "unstructured" (default,
  direct scales and LKJ correlation) or "factor", which uses the sparse
  factor model \\G = \mathrm{diag}(s_G)(\Lambda\Lambda^T +
  \Psi)\mathrm{diag}(s_G)\\ on the correlation scale. The factor model
  gives loadings regularized-horseshoe priors for automatic shrinkage;
  loadings are not exactly zero. The ordinary factor model reports
  `loadings`, `factor_magnitude`, `uniqueness`, `communality` and `s_G`.
  With `residual_structure = "shared_factor"` the same loadings also
  structure E, and the fit reports `loadings`, `factor_magnitude`,
  `factor_heritability`, `genetic_uniqueness`, `residual_uniqueness`,
  `uniqueness` and `communality`. It is available only for multivariate
  fits and is not compatible with `model_file`.

- n_factors:

  Number of factors when `genetic_structure = "factor"`; for
  `residual_structure = "shared_factor"` these are shared by `G` and
  `E`. It is a whole number between 1 and the number of traits minus
  one. If `NULL`, the factor model falls back to
  `min(K - 1, max(2, floor(sqrt(K))), 20)` and warns with the value;
  that fallback depends only on the number of traits, so use
  [`suggest_n_factors`](https://diogro.github.io/stanimal/reference/suggest_n_factors.md)
  to choose it from the data when the factors themselves are meant to be
  interpreted. Ignored otherwise.

- hs_global_scale, hs_slab_scale:

  Positive scales of the regularized-horseshoe prior on factor loadings:
  the global half-Cauchy scale and the slab scale, which limits the
  loading's prior standard deviation rather than bounding its realized
  value. For ordinary factor G, communality is \\q/(1+q)\\, where \\q\\
  is the random sum of squared unnormalized loadings. Its prior
  distribution depends on both scales, the local scales, the number of
  factors and column decay; it is not determined by the slab alone. The
  shared-factor model does not row-normalize its loadings, so its
  communality is a sum of squares on the standardized response scale,
  not a correlation-scale proportion. Report the settings, especially
  where data are weak. These arguments are only used by
  `genetic_structure = "factor"`; explicitly supplying them elsewhere
  warns. See Piironen and Vehtari (2017) on choosing the global scale.

- column_decay:

  Standard-deviation multiplier for successive factor columns, between 0
  and 1. Column `m` is multiplied by `column_decay^(m - 1)`, so values
  below 1 favour decreasing factor magnitudes and make the factor labels
  interpretable. The default 0.8 orders the columns without forcing a
  strict hierarchy; with `column_decay = 1` the columns are exchangeable
  a priori and their labels are free to swap between chains. How strong
  the decay needs to be depends on how many columns have to be told
  apart: with two columns of comparable size 0.8 is often too mild and
  0.5 is what separates them, while a strong decay works against a model
  with many columns, because `0.5^(m - 1)` cuts the prior of the fifth
  column to a sixteenth of the first. This is a deterministic prior on
  the column order, not a learned decay parameter; it is used by both
  factor structures, including the shared-factor model. Explicitly
  supplying it outside `genetic_structure = "factor"` warns that it is
  ignored.

- draws:

  Number of draws from the approximation used for its reported means,
  standard deviations and median absolute deviations. A whole number.

- seed:

  Optional seed for optimization and approximation draws. If omitted, a
  seed is drawn from R's random-number generator and recorded in the
  result.

- correlation_prior:

  Prior on the genetic and residual correlation matrices of the
  unstructured multivariate model. `"lkj"` (default) is
  LKJ(`lkj_prior`), which shrinks every correlation toward zero.
  `"size"` writes each correlation matrix as one dense size factor plus
  an LKJ remainder, \\\Omega = D(ll^T + R)D\\ with \\R \sim\\
  LKJ(`lkj_prior`) and \\D\\ the row normalisation, and estimates where
  the correlations are shrunk to: the size loadings \\l_k = \mu + \sigma
  z_k\\ are pooled around a mean loading with \\\rho = \mu^2/(1+\mu^2)
  \sim\\ Beta(1, 1), \\\sigma \sim\\ half-normal(0, 1) and \\z_k \sim
  N(0, 1)\\. `G` and `E` each get their own. `"shared_size"` is the same
  with one standardized pattern \\z\\ shared by `G` and `E`, so the
  genetic and residual size axes deviate from equal loadings in the same
  direction, each with its own mean loading and spread; it needs
  `residual_structure = "unstructured"`. Neither lowers the rank: the
  factor moves the point the prior shrinks toward, and with \\\rho =
  \sigma = 0\\ the prior is the LKJ one. Use them for traits that share
  a size axis (many measurements of one structure), where LKJ pulls
  large correlations, and with them the heritabilities, toward zero.
  Both need `genetic_structure = "unstructured"`, cannot be combined
  with `model_file`, and use `parametrization = "h2"`: leaving
  `parametrization` unset selects it, and asking for `"variance"` is an
  error. The fits add `rho_G`, `rho_E`, `sigma_size_G`, `sigma_size_E`,
  the size loadings on the correlation scale `size_loadings_G` and
  `size_loadings_E`, and `integration_G` and `integration_E`, the
  relative eigenvalue variance of each correlation matrix (0 for no
  correlation, 1 for a single axis). Sampling is slower than with
  `"lkj"`, about three times in the tests.

- ...:

  additional parameters for stan sampling

- covariance_structure:

  Covariance structure for `lm_phen_cov_laplace`, as in
  [`lm_phen_cov`](https://diogro.github.io/stanimal/reference/lm_phen_cov.md):
  `"unstructured"` (default) or `"factor"`.

## Value

An object of class `stanimal_laplace`, a list with

- `summary` a data frame of posterior means and standard deviations for
  the covariance quantities the model reports, on the original data
  scale. Quantities are labelled with the column names of `Y` (and `X`
  for `beta`), e.g. `G[wing,tail]` and `h2[tail]`, as
  [`fit_draws()`](https://diogro.github.io/stanimal/reference/fit_draws.md)
  names sampler fits; an unnamed `Y` keeps Stan's indices (`G[2,1]`).

- `draws` the draws from the approximation, for inspection only, with
  the same labels.

- `y_names`, `x_names` the column names of `Y` and `X`.

- `elapsed` seconds the fit took.

- `call` the matched call.

It is not a posterior sample and the draw-based functions refuse it. Use
it to check that the sampler and a Laplace approximation find the same
region, to screen many trait sets quickly, or to get a covariance when a
full fit is too slow, and not to report a heritability with an interval.

## Details

The approximation is Gaussian on the unconstrained parameter scale,
centred on the mode, with covariance taken from the curvature there.
Quantities on their reported scales are transformations of those draws.
The approximation is poor wherever the posterior is skewed, and the
diagonal of a covariance matrix is the first place that happens, so a
standard deviation for a variance or a heritability is the least
trustworthy number these functions return. In bivariate animal models
with 120 individuals the means came out systematically below the
sampler's, by up to 0.18 on entries of order one, because the mode sits
below the mean when the posterior is skewed, and the standard deviation
of the smallest variance was about 40% too wide.

Only `cmdstanr` can run the approximation; there is no `rstan`
equivalent, so the functions error when it is not installed rather than
falling back to sampling. Optimization and CmdStan failures and warnings
are surfaced; there is no automatic retry. If the Laplace stage fails,
its CmdStan output is included in the error. The mode is written with
high precision to avoid loss of precision when CmdStan reads a
correlation Cholesky factor back in.

## See also

[`lmm_animal`](https://diogro.github.io/stanimal/reference/lmm_animal.md),
[`lm_phen_cov`](https://diogro.github.io/stanimal/reference/lm_phen_cov.md)
