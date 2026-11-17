# Fit a Bayesian phenotypic covariance matrix

Estimates the phenotypic covariance matrix `P` of one or more traits
after accounting for fixed effects, for data with no relatedness or
genomic information. It is the counterpart of
[`lmm_animal`](https://diogro.github.io/stanimal/reference/lmm_animal.md)
without the genetic component, and uses the same priors, scaling and
outputs.

## Usage

``` r
lm_phen_cov(
  Y,
  X = NULL,
  lkj_prior = 2,
  iter = 3000,
  warmup = floor(iter/2),
  chains = 4,
  cores = chains,
  control = list(adapt_delta = 0.9),
  cmd_stan = TRUE,
  model_file = NULL,
  covariance_structure = c("unstructured", "factor"),
  n_factors = NULL,
  hs_global_scale = 1,
  hs_slab_scale = 1,
  column_decay = 0.8,
  ...
)
```

## Arguments

- Y:

  Numeric response vector or matrix, with individuals in rows. Every
  trait must vary: constant traits cannot be standardized and are an
  error. Missing values (`NA`) are allowed and are not discarded: see
  Details.

- X:

  Optional finite numeric fixed-effects matrix with the same rows as Y
  and full column rank. The default `NULL` fits an intercept-only model:
  one column of ones, named `"(Intercept)"`. A matrix must be able to
  represent a constant (an intercept column, or the 0/1 indicators of
  every level of a factor), so that the mean of the response can be
  reported. Continuous columns are centered and scaled internally;
  returned beta draws use original X units.

- lkj_prior:

  Positive shape of the LKJ prior on the correlation matrix. `1` is
  uniform over matrices, not individual correlations: pairwise marginals
  concentrate increasingly around zero as dimension increases. Values
  above `1` add shrinkage toward zero. Used only by
  `covariance_structure = "unstructured"`; explicitly supplying it for
  the factor model warns that it is ignored.

- iter:

  Number of MCMC iterations, a whole number.

- warmup:

  Number of warmup iterations, a whole number smaller than `iter`.

- chains:

  Number of chains to run, a whole number.

- cores:

  Number of chains to run in parallel (cmdstanr only), a whole number.

- control:

  Named list of sampler controls, passed to whichever engine runs:
  `adapt_delta`, `adapt_engaged`, `max_treedepth`, `step_size` (rstan's
  `stepsize` is accepted), `init_buffer`, `term_buffer`, `window` and
  `metric`. Any other name is an error, so a control is never silently
  ignored.

- cmd_stan:

  Logical, use cmdstanr. If FALSE, uses rstan. Both engines are
  suggested packages: when cmdstanr is not installed the fit falls back
  to rstan with a warning, and when neither is installed the fit errors.

- model_file:

  Optional path to a custom Stan file replacing the default model.

- covariance_structure:

  Covariance structure for the phenotypic model: "unstructured"
  (default, scales times an LKJ correlation matrix) or "factor", which
  uses \\P = \mathrm{diag}(\sigma_P)(\Lambda\Lambda^T +
  \Psi)\mathrm{diag}(\sigma_P)\\ with regularized-horseshoe loadings,
  where \\\sigma_P\\ are the residual standard deviations, estimated
  with the same half-normal prior as the unstructured model.

- n_factors:

  Number of factors when `covariance_structure = "factor"`, a whole
  number between 1 and the number of traits minus one. If `NULL`, the
  model falls back to the same sublinear rule as
  [`lmm_animal()`](https://diogro.github.io/stanimal/reference/lmm_animal.md)
  and warns with the value; that fallback depends only on the number of
  traits, so use
  [`suggest_n_factors`](https://diogro.github.io/stanimal/reference/suggest_n_factors.md)
  to choose it from the data when the factors themselves are meant to be
  interpreted. Ignored otherwise.

- hs_global_scale, hs_slab_scale:

  Positive scales of the regularized-horseshoe prior on the factor
  loadings, used only by the factor model. The slab limits the prior
  standard deviation, not the realized loading. Communality is
  \\q/(1+q)\\, where \\q\\ is the random sum of squared unnormalized
  loadings; its prior depends on both scales, local scales, factor count
  and column decay.

- column_decay:

  Standard-deviation multiplier for successive factor columns, between 0
  and 1. The default 0.8 orders the columns so that the factor labels
  stay meaningful; with `column_decay = 1` the columns are exchangeable
  a priori and their labels can swap between chains. Few factors need a
  stronger decay than many: with two columns of comparable size 0.8 is
  often too mild to separate them, while `0.5^(m - 1)` would shrink the
  fifth column of a many-factor model to a sixteenth of the first. Only
  used by the factor model, and explicitly supplying it otherwise warns
  that it is ignored.

- ...:

  additional parameters for stan sampling.

## Value

With `cmd_stan = TRUE`, a `CmdStanMCMC` object from cmdstanr, usable
with posterior and bayesplot (see
[`lmm_animal`](https://diogro.github.io/stanimal/reference/lmm_animal.md));
use
[`rename_beta`](https://diogro.github.io/stanimal/reference/rename_beta.md)
for `beta` draws named after the columns of `X`, and
[`fit_draws`](https://diogro.github.io/stanimal/reference/fit_draws.md)
or
[`draws_mean`](https://diogro.github.io/stanimal/reference/draws_mean.md)
to get `P`, `corrP` and `sigma_P` as matrices. With `cmd_stan = FALSE`,
a `stanfit` object.

The unstructured fit reports `P` (original trait scale), `corrP`,
`sigma_P` and the fixed effects `beta[trait, j]`. `sigma_P` are the
residual standard deviations, so `P`, `corrP` and `sigma_P` describe the
covariance of the residuals after the fixed effects, not the marginal
covariance of `Y`.

Factor fits add quantities that describe the decomposition.
`communality` and `uniqueness` do not depend on how the factors are
oriented, so they are identified by the data. `loadings` and
`factor_magnitude` do: two solutions with the same `P` can differ by a
rotation, a sign flip or a permutation of the columns, and nothing fixes
that choice except the prior. A poor `Rhat` on those two is therefore
not evidence that the fit failed, and their posterior means describe
whichever orientations the chains visited. Judge convergence on `P`,
`corrP`, `sigma_P`, `communality` and `uniqueness`, and pin the
orientation down by lowering `column_decay`, which orders the columns,
or by asking for fewer factors.

- `loadings`: the trait by factor loading matrix \\\Lambda\\, on the
  correlation scale.

- `communality`: for trait \\k\\, \\c_k = \sum_m \lambda\_{km}^2\\, the
  part of the trait's variance explained by the common factors. On the
  correlation scale \\c_k + u_k = 1\\; multiplying by `sigma_P[k]^2`
  gives the common phenotypic variance of trait \\k\\.

- `uniqueness`: for trait \\k\\, \\u_k = 1 - c_k\\, the trait-specific
  variance left after the factors. Multiplying by `sigma_P[k]^2` gives
  the unique phenotypic variance.

- `factor_magnitude`: for factor \\m\\, \\\sum_k \lambda\_{km}^2\\, the
  variance the factor contributes across traits. It is a factor-level
  summary, not a per-trait quantity, and its values only compare across
  columns once the orientation is fixed.

These are available from `fit$summary()` and from
[`fit_draws`](https://diogro.github.io/stanimal/reference/fit_draws.md):
`loadings` is a trait by factor by draw array, `factor_magnitude` is a
factor by draw matrix, and the uniqueness and communality vectors are
trait by draw matrices.

## Details

The unstructured model estimates `P` as scales times an LKJ-distributed
correlation matrix, with half-normal priors on the scales of the
standardized responses. With `covariance_structure = "factor"`, `P` is
instead \\\mathrm{diag}(\sigma_P)(\Lambda\Lambda^T +
\Psi)\mathrm{diag}(\sigma_P)\\, with regularized-horseshoe loadings and
deterministic column decay. The factor fit reports `loadings`,
`factor_magnitude`, `uniqueness` and `communality` in addition to `P`,
`corrP`, `sigma_P` and the fixed effects. Returned quantities are on the
original response scale.

Incomplete records are used rather than dropped. Without relatedness the
records are independent, so a row that misses some traits is evaluated
on the traits it does have: its marginal density is the multivariate
normal on that subvector, with the matching submatrix of `P`. Rows are
grouped by which traits they observe, each group gets its own Cholesky
factor, and the missing values are integrated out analytically, so they
cost no parameters. A record with nothing observed carries no likelihood
and is an error instead. The count of groups grows with the number of
distinct patterns, at most one per row.

The pointwise `log_lik` is the density the model block evaluates, on the
standardized response; the raw-scale density differs from it by a
constant, which does not affect comparisons between fits to the same
data.

The likelihood is computed by whitening the residuals with the Cholesky
factor of `P` in a single triangular solve, at `O(N K^2 + K^3)` per
evaluation rather than a per-individual loop, so the model scales to
dozens of traits.

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
n <- 60
X <- cbind(Intercept = 1, sex = rbinom(n, 1, 0.5))
P <- random_covariance(3, min_var = 0.8, max_var = 1.2,
                       trait_names = c("t1", "t2", "t3"))
B <- rbind(c(10, 9, 11), c(1, -0.5, 0.3))
Y <- X %*% B + matrix(rnorm(n * 3), n) %*% chol(P)
colnames(Y) <- c("t1", "t2", "t3")

fit <- lm_phen_cov(Y, X, iter = 1000, warmup = 500,
                   chains = 2, cores = 2, refresh = 0)
fit$summary(c("P", "corrP", "sigma_P"))
round(draws_mean(fit, "P"), 2)

# Sparse factor P
factor_fit <- lm_phen_cov(Y, X, covariance_structure = "factor",
                          n_factors = 2, iter = 1000, warmup = 500,
                          chains = 2, cores = 2, refresh = 0)
factor_fit$summary(c("P", "loadings", "factor_magnitude"))
dim(fit_draws(factor_fit)$loadings)
} # }
```
