# Fit linear mixed animal model

Fits a univariate or multivariate linear mixed animal model using a
known relatedness matrix

## Usage

``` r
lmm_animal(
  Y,
  X = NULL,
  A,
  lkj_prior = 2,
  beta_prior = c(1, 1),
  parametrization = c("variance", "h2"),
  iter = 3000,
  warmup = floor(iter/2),
  chains = 4,
  cores = chains,
  model_file = NULL,
  control = NULL,
  cmd_stan = TRUE,
  blend_A = FALSE,
  residual_structure = c("unstructured", "diagonal", "shared_factor"),
  genetic_structure = c("unstructured", "factor"),
  n_factors = NULL,
  hs_global_scale = 1,
  hs_slab_scale = 1,
  column_decay = 0.8,
  eigen_A = NULL,
  correlation_prior = c("lkj", "size", "shared_size"),
  ...
)
```

## Arguments

- Y:

  response variable, univariate or multivariate. Every trait must vary:
  constant traits cannot be standardized and are an error. Missing
  values (`NA`) are allowed and are not discarded: see Details.

- X:

  Optional finite numeric fixed-effects matrix with the same rows as Y
  and full column rank. The default `NULL` fits an intercept-only model:
  one column of ones, named `"(Intercept)"`. A matrix must be able to
  represent a constant (an intercept column, or the 0/1 indicators of
  every level of a factor), so that the mean of the response can be
  reported. Continuous columns are centered and scaled internally;
  returned beta draws use original X units.

- A:

  Finite symmetric positive-definite relationship matrix matching the
  rows of Y. A matrix that is only positive semidefinite (duplicated
  individuals, for example) is rejected; set `blend_A = TRUE` to
  regularize it. When `A` and `Y` both have row names, they are used to
  put `A` in the order of `Y`; names that disagree are an error. Without
  names the rows are matched by position.

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

- iter:

  Number of MCMC iterations, a whole number.

- warmup:

  Number of warmup iterations, a whole number smaller than `iter`.

- chains:

  Number of chains to run, a whole number.

- cores:

  Number of chains to run in parallel (cmdstanr only), a whole number.

- model_file:

  Optional path to a custom Stan file replacing the default model.

- control:

  Named list of sampler controls, passed to whichever engine runs:
  `adapt_delta`, `adapt_engaged`, `max_treedepth`, `step_size` (rstan's
  `stepsize` is accepted), `init_buffer`, `term_buffer`, `window` and
  `metric`. Any other name is an error, so a control is never silently
  ignored. The default `NULL` uses `list(adapt_delta = 0.9)` for
  `genetic_structure = "factor"` with an unstructured or diagonal
  residual, and `list(adapt_delta = 0.99)` otherwise: unstructured `G`
  on small samples and shared-factor fits diverge at 0.9. A list you
  supply replaces the default entirely.

- cmd_stan:

  Logical, use cmdstanr. If FALSE, uses rstan. Both engines are
  suggested packages: when cmdstanr is not installed the fit falls back
  to rstan with a warning, and when neither is installed the fit errors.

- blend_A:

  If TRUE, replace A with 0.99 \* A + 0.01 \* I before fitting, making
  singular A positive definite.

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

- eigen_A:

  Optional eigendecomposition of `A`, a list with elements `Q` (the \\N
  \times N\\ orthonormal eigenvectors, as columns) and `lambda` (the `N`
  matching eigenvalues), for example
  `e <- eigen(A, symmetric = TRUE); list(Q = e$vectors, lambda = e$values)`.
  The models never decompose `A` themselves: without `eigen_A` it is
  decomposed once per fit in R. Supplying it skips that step, so a
  caller that fits many models on the same relationship matrix (one per
  trait subset, say) pays the `O(N^3)` decomposition once. It must be
  the decomposition of `A` in the row order of `A` as passed; rows are
  reordered along with `A` when `A` and `Y` are matched by name, and
  `blend_A` is applied to the eigenvalues. Eigenvalues and eigenvectors
  need only be consistently paired, not in any particular order, though
  the order sets the order of the rotated records that `log_lik` refers
  to. The pair is checked against `A` (reconstruction and
  orthonormality, on two probe vectors) before fitting. Works with every
  structure and parametrization.

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

## Value

With `cmd_stan = TRUE`, a `CmdStanMCMC` object from cmdstanr: use
`fit$draws()` (a `posterior` draws array) and `fit$summary()`, and pass
these, or the fit itself, to bayesplot (e.g.
`bayesplot::mcmc_trace(fit$draws("h2"))`,
`bayesplot::nuts_params(fit)`). Fixed effects are `beta[j]` (univariate)
or `beta[trait,j]` for column `j` of `X`; use
[`rename_beta`](https://diogro.github.io/stanimal/reference/rename_beta.md)
to get draws named after the columns of `X`. Use
[`fit_draws`](https://diogro.github.io/stanimal/reference/fit_draws.md)
to reshape the draws into trait by trait by draw arrays, and
[`draws_mean`](https://diogro.github.io/stanimal/reference/draws_mean.md)
or
[`draws_quantile`](https://diogro.github.io/stanimal/reference/draws_mean.md)
for posterior summaries that are already matrices. Factor fits also
expose `loadings`, `factor_magnitude`, `uniqueness` and `communality`;
shared-factor fits add `factor_heritability`, `genetic_uniqueness` and
`residual_uniqueness`. With `cmd_stan = FALSE`, a `stanfit` object whose
`beta` draws are named after the columns of `X`. Multivariate fits
report `G`, `E`, `P = G + E`, `corrG`, `sigma_G`, `sigma_E` and
`h2 = diag(G) / diag(P)`, plus the fixed effects `beta[trait, j]`.

Every fit reports `log_lik`, one term per rotated response record. For
complete data the terms sum to the joint log density, but leave-one-out
over rotated records is not generally leave-one-individual-out.
Comparisons need the same response and predictive units, including the
same relationship matrix and eigenbasis. With missing responses the
terms use the completed data rather than the marginalized observed-data
likelihood, so
[`fit_loo`](https://diogro.github.io/stanimal/reference/fit_loo.md)
refuses these fits.

Factor fits add quantities that describe the factor decomposition. They
come in two groups. `communality`, `uniqueness` and
`factor_heritability` do not depend on how the factors are oriented, so
the data identify them and their `Rhat` and effective sample size are
read like any other parameter. `loadings` and `factor_magnitude` do
depend on it: two solutions with the same `G` can differ by a rotation,
a sign flip or a permutation of the columns, and nothing fixes that
choice except the prior. A poor `Rhat` on those two is therefore not
evidence that the fit failed, and their posterior means describe
whichever orientations the chains visited. Judge convergence on `G`,
`E`, `P`, `corrG`, `communality` and `uniqueness`, and pin the
orientation down by lowering `column_decay`, which orders the columns,
or by asking for fewer factors.

- `loadings`: the trait by factor loading matrix \\\Lambda\\. For the
  ordinary factor model these are on the correlation scale; for
  shared-factor fits they are on the standardized phenotype scale.

- `communality`: for trait \\k\\, \\c_k = \sum_m \lambda\_{km}^2\\. This
  is the part of the trait's variance explained by the common factors.
  In the ordinary factor model it is on the correlation scale, so
  \\c_k + u_k = 1\\; multiplying by \\s\_{Gk}^2\\ gives the common
  genetic variance of trait \\k\\. In shared-factor fits it is the
  common phenotypic variance, and \\c_k + u_k\\ equals the trait's
  standardized phenotypic variance.

- `uniqueness`: for trait \\k\\, the trait-specific variance left after
  the factors. In the ordinary factor model \\u_k = 1 - c_k\\ is the
  unique genetic variance on the correlation scale; multiply by
  \\s\_{Gk}^2\\ for the unique genetic variance. In shared-factor fits
  \\u_k = \psi\_{Gk} + \psi\_{Ek}\\, with `genetic_uniqueness`
  \\\psi\_{Gk}\\ and `residual_uniqueness` \\\psi\_{Ek}\\ giving the two
  parts.

- `factor_magnitude`: for factor \\m\\, \\\sum_k \lambda\_{km}^2\\, the
  variance the factor contributes across traits. It is a factor-level
  summary, not a per-trait quantity. Ranking the columns only means
  something once the orientation is fixed, by `column_decay` or by an
  `n_factors` small enough to leave the factors identified.

- `factor_heritability`: shared-factor fits only. For factor \\m\\,
  \\h^2_m\\ is the fraction of that factor's variance that is genetic:
  the genetic part is multiplied by \\h^2_m\\ and the residual part by
  \\1-h^2_m\\.

- `s_G`: the genetic standard deviation of each trait on the
  standardized response scale, which is the diagonal scale that moves
  the correlation-scale decomposition to the original scale; multiplying
  it by the response standard deviation gives the genetic standard
  deviations `sigma_G` on the data scale. It is reported by every
  multivariate fit except the shared-factor model, which scales its
  factors by `sqrt_psi_G` and `sqrt_psi_E` instead. In the ordinary
  factor model it is derived from `h2` and `s_P`.

These are available from `fit$summary()` and from
[`fit_draws`](https://diogro.github.io/stanimal/reference/fit_draws.md).
The latter returns `loadings` as a trait by factor by draw array,
`factor_heritability` and `factor_magnitude` as factor by draw matrices,
and the uniqueness and communality vectors as trait by draw matrices.

## Details

By default the models estimate the additive genetic and residual
standard deviations (`sigma_G`, `sigma_E`) directly, with half-normal
priors on the standardized response scale, and report `h2` as a derived
quantity. With `parametrization = "h2"`, and always in the ordinary
factor model, the two scales are written as \\s_G = \sqrt{h^2}s_P\\ and
\\s_E = \sqrt{1-h^2}s_P\\ instead, so `beta_prior` is a prior on the
trait heritabilities and \\s_P\\ is the phenotypic scale. The
multivariate model estimates the additive genetic covariance matrix `G`
and the residual covariance `E` directly, each as scales times an
LKJ-distributed correlation matrix (`lkj_prior`), with half-normal
priors on the scales of the standardized responses. With
`genetic_structure = "factor"`, `G` is instead
\\\mathrm{diag}(s_G)(\Lambda\Lambda^T + \Psi)\mathrm{diag}(s_G)\\, where
the loadings \\\Lambda\\ have regularized-horseshoe priors and \\\Psi\\
is diagonal on the correlation scale; the default global shrinkage is
`hs_global_scale = 1`. With `residual_structure = "shared_factor"`, the
same loadings enter both covariances: \\G =
\Lambda\mathrm{diag}(h^2\_\text{factor})\Lambda^T +
\mathrm{diag}(\psi_G)\\ and \\E = \Lambda\mathrm{diag}(1 -
h^2\_\text{factor})\Lambda^T + \mathrm{diag}(\psi_E)\\, with the beta
prior on the factor heritabilities controlled by `beta_prior`. The
optional `column_decay` multiplies the prior standard deviation of
factor column \\m\\ by `column_decay^(m - 1)`, which favours factors in
decreasing order without adding a learned decay parameter. Breeding
values are marginalized using the eigendecomposition of `A`: the
response and the design are rotated by its eigenvectors and the
likelihood is computed on the rotated records, with the eigenvalues
entering through the per-record covariance. No per-individual random
effects are needed, which is why a singular `A` can be fitted once it is
regularized with `blend_A = TRUE`. Returned quantities (original
response scale) are `G`, `E`, `P = G + E`, `corrG`, `sigma_G`, `sigma_E`
and `h2 = diag(G) / diag(P)`.

With `parametrization = "h2"`, and always in the ordinary factor model,
the scales are written as a phenotypic scale and a heritability rather
than as two standard deviations, because the response is standardized
trait by trait and the two scales are tied to its variance. That makes
\\s_G^2 + s_E^2\\ the phenotypic variance of the trait and fixes it near
1, so half-normal priors on \\s_G\\ and on \\s_E\\ are priors on a
circle. The density they induce on \\h^2\\ is proportional to
\\1/\sqrt{h^2(1-h^2)}\\, an arcsine prior that puts about 14 percent of
its mass below 0.05 and 14 percent above 0.95, which makes a moderate
heritability the least likely answer a priori. With few individuals the
likelihood is weak enough for that to show through, and a fit can then
sit at either end of the interval. Writing the scales as \\s_G =
\sqrt{h^2}s_P\\ and \\s_E = \sqrt{1-h^2}s_P\\ puts `beta_prior` on
\\h^2\\ directly, so the prior applied is the one that was asked for,
and \\s_P\\ is the well-determined quantity that is left.

The per-individual likelihoods are whitened by the residual Cholesky
factor and rotated into a single eigendecomposition of
`E^{-1/2} G E^{-1/2}`, so each evaluation costs `O(N K^2 + K^3)` rather
than `O(N K^3)`. This lets the multivariate model scale to dozens of
traits without a Cholesky factorization per individual. When the
eigenvalues of that matrix are nearly coincident (for example at
initialization, or when `G` is roughly proportional to `E`) the model
automatically falls back to a per-individual Cholesky evaluation, which
has a well-behaved gradient.

The cost that grows with the number of individuals is the dense
relationship matrix: `A` is passed in full and each chain
eigendecomposes it before sampling (`O(N^3)` time and `O(N^2)` memory
per chain). That cost is paid once, so with many traits sampling
dominates and the decomposition is a small share of the fit, while with
few traits and short runs it can be most of it. The models are intended
for many traits on hundreds to a few thousand individuals. For cohorts
of tens of thousands, where that matrix is the bottleneck, an engine
that uses the sparse inverse of `A` (such as ASReml) is the better tool.

Incomplete records are used rather than dropped. A missing value is not
discarded, because with relatedness it carries information: the value is
informed by every relative's observed traits. The gaps are filled by
parameters on the standardized response scale, drawn with the rest of
the model and reported as `Y_mis`, and a flat prior on them makes the
posterior over the other parameters the exact marginal over the observed
data. Each gap is one parameter, so the cost grows with the number of
gaps; a record with nothing observed is allowed, since its relatives
still inform it. Dropping incomplete rows instead is the alternative,
and it throws away the rows that are complete for the other traits.

Two configurations are worth avoiding when many values are missing. A
`residual_structure = "diagonal"` model with a small sample is weakly
identified once gaps loosen the residual structure, and can mix badly:
prefer `"unstructured"`, which is the default and what the missing-data
machinery is tested against. For fits with gaps, `log_lik` uses the
completed response including imputation parameters, not the marginalized
observed-data likelihood;
[`fit_loo`](https://diogro.github.io/stanimal/reference/fit_loo.md)
refuses LOO and WAIC for these fits. Raw `log_lik` draws remain
available, but are not an observed-data predictive criterion.

The LKJ shape `lkj_prior` (default 2) deliberately shrinks correlations
toward zero. With many traits, an unstructured correlation matrix fit
under a flat prior tends to overestimate the magnitude of correlations,
especially genetic correlations between traits with limited information,
so a shape above 1 acts as a mild regularizer; increase it for stronger
shrinkage. Note that fitting many traits remains expensive: the
likelihood is cheap, but an unstructured `G` for `K = 40` traits has on
the order of 1700 parameters and NUTS needs many leapfrog steps, so a
single fit can take an hour or more. Fewer iterations, a larger
`lkj_prior`, or `residual_structure = "diagonal"` all make fits faster.

Continuous predictors are centered and divided by their sample standard
deviation inside Stan. Constant and 0/1 indicator columns are not
centered or scaled. Standardized slopes have `normal(0, 1)` priors; the
intercept contribution has a wider `normal(0, 2.5)` prior. Returned
`beta` coefficients are backtransformed to original predictor and
response units, including restoring the response mean to the intercept;
`beta_scaled` draws are also available.

## Author

Diogo Melo

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
n <- 40
ped <- data.frame(
  id = 1:n,
  sire = c(rep(NA, 20), sample(1:10, 20, replace = TRUE)),
  dam = c(rep(NA, 20), sample(11:20, 20, replace = TRUE))
)
A <- makeA(ped)
X <- cbind(Intercept = 1, sex = rbinom(n, 1, 0.5))
G <- random_covariance(2, min_var = 0.8, max_var = 1.2,
                       trait_names = c("t1", "t2"))
E <- random_covariance(2, min_var = 0.8, max_var = 1.2,
                       trait_names = c("t1", "t2"))
Y <- X %*% rbind(c(10, 9), c(1, -0.5)) +
  t(chol(A)) %*% matrix(rnorm(n * 2), n) %*% chol(G) +
  matrix(rnorm(n * 2), n) %*% chol(E)
colnames(Y) <- c("t1", "t2")

fit <- lmm_animal(Y, X, A, iter = 1000, warmup = 500,
                  chains = 2, cores = 2, refresh = 0)
fit$summary(c("G", "E", "h2"))
round(draws_mean(fit, "G"), 2)

# Univariate model with the h2 parametrization
fit_h2 <- lmm_animal(Y[, 1], X, A, parametrization = "h2",
                     iter = 1000, warmup = 500,
                     chains = 2, cores = 2)

# Sparse factor G
factor_fit <- lmm_animal(Y, X, A, genetic_structure = "factor",
                         n_factors = 1, iter = 1000, warmup = 500,
                         chains = 2, cores = 2)
factor_fit$summary(c("loadings", "factor_magnitude"))
dim(fit_draws(factor_fit)$loadings)

# Shared factors for G and E
shared_fit <- lmm_animal(Y, X, A, genetic_structure = "factor",
                         residual_structure = "shared_factor",
                         n_factors = 1, iter = 1000, warmup = 500,
                         chains = 2, cores = 2)
shared_fit$summary(c("factor_heritability", "genetic_uniqueness"))
} # }
```
