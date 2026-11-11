# Using stanimal — guide for agents

`stanimal` fits Bayesian quantitative-genetic models in Stan (CmdStan via cmdstanr,
with an rstan fallback). Covariances, standard deviations, heritabilities and
fixed effects are returned on the **original data scale**, so they can be compared
with the trait values directly. In ordinary factor models, `loadings` and
`communality`/`uniqueness` are on the correlation scale; in shared-factor
models they are on the standardized response scale. `s_G` is on the
standardized response scale.

This file is the interface guide: what to call, what comes back, and how to make
the prior choices that the function signatures do not make for you. It is written
for anyone *using* the package. Worked examples live in `README.md` and
`examples/`; if you are changing the package itself, read `DEVELOPMENT.md`.

## Quick decision table

| Goal | Call |
|---|---|
| One trait, heritability | `lmm_animal(y, X, A)` |
| Several traits, covariance structure | `lmm_animal(Y, X, A)` |
| Many traits, `G` needs to be sparse/low-rank | `lmm_animal(Y, X, A, genetic_structure = "factor")` |
| `G` and `E` should share the same latent factors | add `residual_structure = "shared_factor"` (needs `"factor"`) |
| No pedigree, only the phenotypic covariance | `lm_phen_cov(Y, X)` |
| Several traits, some values missing | pass the `NA`s straight in — see below |
| The completed data, for follow-up analyses | `imputed_data(fit)` |
| Phenotypic covariance of many traits | `lm_phen_cov(Y, X, covariance_structure = "factor")` |
| Build `Y`/`X` from a formula, with `A` | `animal_design(formula, data, pedigree = ...)` |
| Pedigree with related founders/metafounders | `makeA(pedigree, Gamma = , groups = )` |
| Markers instead of a pedigree | `genomic_relationship(M)` → pass as `A` |
| How many factors? | `suggest_n_factors(Y, X)` + `scree_eigenvalues(Y, X)` |
| Which of two models fits better? | `loo::loo_compare(fit_loo(a), fit_loo(b))` (a size-shape fit takes `fit_loo(b, transform = ss)`) |
| Quick approximate fit (no draws-based tools) | `lmm_animal_laplace(Y, X, A)` / `lm_phen_cov_laplace(Y)` |
| Does the model reproduce the data? | `posterior_predictive(fit)` |
| Each individual's breeding value | `breeding_values(fit)` |
| Many correlated traits on one size gradient | `size_shape()` then `back_transform()` — but only around an **unstructured** model, see the anti-patterns |

`A` is the third positional argument of `lmm_animal()`, so `lmm_animal(Y, X, A)`.

## The interface

```r
lmm_animal(Y, X = NULL, A, lkj_prior = 2, beta_prior = c(1, 1),
           parametrization = c("variance", "h2"),
           iter = 3000, warmup = floor(iter/2), chains = 4, cores = chains,
           model_file = NULL, control = NULL,
           cmd_stan = TRUE, blend_A = FALSE,
           residual_structure = c("unstructured", "diagonal", "shared_factor"),
           genetic_structure = c("unstructured", "factor"), n_factors = NULL,
           hs_global_scale = 1, hs_slab_scale = 1, column_decay = 0.8,
           eigen_A = NULL, ...)

lm_phen_cov(Y, X = NULL, lkj_prior = 2, iter = 3000, warmup = floor(iter/2),
            chains = 4, cores = chains, control = list(adapt_delta = 0.9),
            cmd_stan = TRUE, model_file = NULL,
            covariance_structure = c("unstructured", "factor"), n_factors = NULL,
            hs_global_scale = 1, hs_slab_scale = 1, column_decay = 0.8, ...)

lmm_animal_gibbs(Y, X = NULL, A, iter = 1000, warmup = floor(iter/2), thin = 1)
```

`Y` is a numeric matrix (individuals × traits) or a vector; `X` defaults to an
intercept; `A` is a symmetric positive-definite relationship matrix whose rows
are matched to `Y` by name when both have names. `A` is decomposed once per fit
(`eigen(A)`, the `O(N^3)` step); when fitting many models on the same `A`, pass
`eigen_A = list(Q = e$vectors, lambda = e$values)` with
`e <- eigen(A, symmetric = TRUE)` to do it once.

`control = NULL` means `adapt_delta = 0.9` for `genetic_structure = "factor"`
with an unstructured or diagonal residual, and 0.99 for everything else
(univariate, unstructured `G`, shared factor), which diverge at 0.9. A `control` list you pass replaces that default entirely, so
include `adapt_delta` in it if you set any other control.

### The relationship matrix enters through its eigendecomposition (since 0.8.0)

The animal models never decompose `A` in Stan any more: `lmm_animal()` computes
`eigen(A)` once per fit in R (in Stan's ascending order) and passes `Q` and
`lambda` as data. What that means in practice:

- **Fitting many models on one `A`** — trait subsets, bootstrap replicates,
  prior sensitivity runs — decompose it once and pass `eigen_A` to every fit.
  At N = 4000 that removes about a minute of setup from each fit.
- **Speed follows R's BLAS.** With the reference LAPACK a large single fit is
  slightly slower than before (~6% of wall clock at N = 4000); with OpenBLAS,
  MKL or Accelerate it is faster. `eigen_A` also accepts a decomposition from
  any other solver (e.g. `torch::linalg_eigh()` on a GPU).
- **`eigen_A` is checked, not trusted**: a mismatched or non-orthonormal pair
  is an error before anything is fitted. Pass the decomposition of `A` in the
  row order you pass `A`; reordering by name and `blend_A` are applied for you.
- **The model did not change.** Same likelihood, same priors, and in the tests
  draw-for-draw identical fits to 0.7.x. With repeated eigenvalues in `A`
  (common for pedigrees) the basis within an eigenspace may differ from
  0.7.x's, so a fixed seed can give a different trajectory to the same
  posterior. The model's `log_lik` terms still refer to records in the
  (ascending) eigenbasis of `A`; `fit_loo()` no longer uses them by default
  (see Model comparison).
- A custom `model_file` written for 0.7.x that declares `A` keeps working: it
  receives `A` as well as `Q` and `lambda`. `...` is forwarded to the
sampler, so `seed`, `refresh`, `init` and `show_messages` work there.

Both model functions return a `CmdStanMCMC` (or an `stanfit` with
`cmd_stan = FALSE`). A fit records the configuration it ran with, and the two
functions record different things:

- `lmm_animal()`: `genetic_structure`, `residual_structure`, `parametrization`,
  `n_factors` (factor fits), `lkj_prior`, `beta_prior`.
- `lm_phen_cov()`: `covariance_structure`, `n_factors` (factor fits),
  `lkj_prior`, `column_decay`.
- Both: `x_names`, `y_names`, and `custom_model` (TRUE when `model_file` was
  used). A custom file can use any prior, so nothing about the prior can be
  recovered from those attributes.

### What each configuration reports

| fit | quantities |
|---|---|
| univariate | `sigma_G`, `sigma_E`, `h2`, `beta`, `log_lik` |
| multivariate, unstructured | `G`, `E`, `P`, `corrG`, `sigma_G`, `sigma_E`, `h2`, `beta`, `log_lik` |
| factor `G` | + `loadings`, `uniqueness`, `communality`, `factor_magnitude`, `s_G` |
| shared factor | + `loadings`, `uniqueness`, `communality`, `factor_magnitude`, `factor_heritability`, `genetic_uniqueness`, `residual_uniqueness` |
| `lm_phen_cov()` | `P`, `corrP`, `sigma_P`, `beta`, `log_lik` (+ `loadings`, `uniqueness`, `communality`, `factor_magnitude` for `"factor"`) |

`h2` exists only for animal fits. Every model reports `log_lik`, one term per
record. For complete animal responses the terms sum to the joint log density,
but refer to records in the eigenbasis of `A`, not generally to individuals.
Incomplete animal fits instead use completed-data densities. `fit_loo()` does
not use these animal terms by default; it computes one term per individual
after the fit (see Model comparison).

### Reading the draws

```r
fit$summary(c("G", "h2"))          # raw variables, tidy data frame
fit_draws(fit)                     # list of arrays: traits first, draws last
draws_mean(fit, "G")               # posterior mean as a named matrix
draws_quantile(fit, "corrG")       # 2.5% and 97.5% as a 3-d array
rename_beta(fit)                   # beta named after the columns of X
```

For a fit of size/shape scores, `to_draws_list(back_transform(fit, transform))`
returns one posterior draws matrix per back-transformed quantity; `transform`
is the corresponding `size_shape()` object.

`draws_mean()`/`draws_quantile()` take `parameter = NULL` for everything, and a
`traits` argument to relabel. Factor extras keep their natural shape:
`loadings` is **trait × factor × draw**, and `factor_magnitude` and
`factor_heritability` are **factor × draw**.

## Incomplete trait records

Put the `NA`s in `Y` and fit. Nothing is dropped: a row survives if it has at
least one observed trait, and the missing values are handled by the likelihood.
Every trait needs at least two observed values, since it has to be standardized.

The two engines do it differently, for a reason that matters when reading output.

- **`lm_phen_cov()` — per pattern.** Without relatedness the records are
  independent, so a row's missing values can be integrated out analytically: its
  density is the normal on the observed subvector, with the matching submatrix of
  `P`. Rows are grouped by which traits they observe (a *pattern*), one Cholesky
  per pattern, no extra parameters. A row with nothing observed carries no
  likelihood and is an error.
- **`lmm_animal()` — imputation parameters.** With relatedness a gap is informed
  by every relative's observed traits, so the marginal does not factorize over
  rows and there is no cheap per-pattern factor. The gaps become parameters
  (`Y_mis`, on the standardized response scale) with a flat prior, which makes
  the posterior over `G`, `E` and `beta` the exact marginal over the observed
  data.

Either way the completed data come out through `imputed_data(fit)`, not from the
raw draws: `fit$draws("Y_mis")` is on the standardized scale and carries no key
to which individual and trait each value belongs to. `imputed_data()` returns
individuals by traits by draws with the gaps filled on the response scale, and
draws the phenotypic model's gaps from their conditional distribution, since that
model has no gap parameters to report. These are multiple imputations: do not
collapse them to one filled-in matrix and treat it as data.

Practical consequences:

- `lmm_animal()` has one parameter per gap. That is fine at a few hundred gaps
  and noticeable at tens of thousands; `lm_phen_cov()` adds none.
- A fully missing record is allowed for `lmm_animal()` (its relatives inform it)
  but not for `lm_phen_cov()`.
- Drop rows with missing *predictors*: `animal_design()` still does that, since a
  row with an unknown fixed effect cannot be placed in the design. Missing
  responses are kept.
- `residual_structure = "diagonal"` with a small sample and many gaps is weakly
  identified and can mix badly. Use `"unstructured"` (the default) when gaps are
  present.
- The gap parameters are latent, so a high `Rhat` on `Y_mis` says little on its
  own; judge convergence on `G`, `E`, `P` and `h2` as usual.

Before the traits reach a model they may pass through something else, and those
take gaps too:

- `size_shape()` completes the response first, because neither size (a mean of
  the traits) nor shape (a difference between them) is defined for a row with a
  gap. The completion is a maximum-likelihood multivariate normal (EM), and a
  gap becomes its conditional mean given the traits that individual has. The
  object records `gaps` and `scores_sd`, the latter being the standard deviation
  of each completed score: it is `0` for rows that were fully measured, and it is
  the reason a score should not be compared across individuals without it.
- The default requests `ncol(Y) - 1` shape axes. If the centered completed shape
  data have insufficient numerical rank, set `n_shape` to at most the reported
  rank; no more than `min(nrow(Y) - 1, ncol(Y) - 1)` axes can vary.
- EM errors if it does not converge or its covariance is singular.
  `animal_design()` still returns the design but warns if the PC diagnostic
  cannot be computed.
- This is *point* completion, unlike `imputed_data()`. The scores are then
  passed to a model that treats them as data, so the uncertainty in
  `scores_sd` does not reach that model's posterior. Say so if asked why the
  intervals from a fit on completed scores look narrow.
- `pc1_variance_share()` and `plot_G_ellipses()` take gaps as well: the share
  comes from the same EM covariance, and a panel of the ellipse plot draws only
  the individuals that have both of its traits. `animal_design()` keeps a
  response with gaps and warns about the dominant component from the incomplete
  data rather than from the complete rows only.

## Model comparison and predictive checks

`fit_loo()` leaves out one individual at a time, for every built-in model and
with or without missing traits.

```r
fit_loo(fit)                                   # PSIS-LOO
loo::loo(fit)                                  # same thing, via the method
loo::loo_compare(fit_loo(a), fit_loo(b))       # differences in ELPD
posterior_predictive(fit)                      # replicates of the response
```

- `fit_loo()` hands pointwise terms to `loo::loo`, with the chain structure
  kept so the effective-sample-size correction works. `criterion = "waic"` gives
  WAIC instead. It needs the `loo` package (a Suggests).
- **Animal fits: one term per individual (`unit = "individual"`, the default).**
  Relatives share breeding values, so the model does not factorize over
  animals. The term for animal `i` is the density of its observed traits given
  every other animal, computed after the fit from the draws of `G`, `E`,
  `beta` and, with gaps, `Y_mis` (non-factorized LOO, Bürkner, Gabry & Vehtari
  2021). The animal's own gaps are integrated out; the other animals' gaps are
  conditioned on at their sampled values, which keeps the importance sampling
  exact. An animal with nothing observed has no term. `attr(result, "units")`
  names the animal of each term, so a high Pareto k can be traced.
- Cost: one `eigen(A)` (pass `eigen_A` to skip it) plus two `N x N` by
  `N x K` products per draw. About 45 s for 505 animals, 12 traits and 8,000
  draws; gaps add a little.
- `unit = "record"` uses the model's own `log_lik` instead: one term per record
  rotated by the eigenvectors of `A`, which is not one animal. It is refused
  for responses with gaps, whose `log_lik` uses imputed values.
- Phenotypic fits already have one term per individual, with gaps integrated
  out, and use them directly.
- Compare fits to the same response values, observations and missingness, with
  the same `unit`. Individual-level terms assume the built-in likelihood and
  are refused for a custom `model_file`. Bare matrices and older fits without
  response metadata fall back to their `log_lik` with a warning. Passing
  extracted draws directly to `loo::loo()` or `loo::waic()` bypasses the
  checks.
- The terms are the density of the *standardized* response, so the absolute ELPD
  carries a constant that depends on the data. It cancels between fits to the
  same response, which is what `loo_compare()` uses.
- **Raw traits against size and shape scores.** A fit to `size_shape()` scores
  is a fit to a linear map of the traits, so its terms can be moved to the
  trait scale: `loo_compare(fit_loo(raw), fit_loo(ss_fit, transform = ss))`.
  The shift is one constant (the Jacobian of the map and the change of
  standardization), exact for a complete response with every shape axis kept,
  and refused otherwise. Without `transform` the two are densities of
  different responses and must not be compared.
- Check `fit_loo(fit)$diagnostics$pareto_k` before trusting a comparison. With a
  few dozen individuals, k above 0.7 is common and `loo` says so itself; the fix
  is a longer run or reading the comparison as indicative.
- `posterior_predictive()` is the individual-level check: it simulates the
  response from each of `ndraws` posterior draws using the reported `G`, `E`
  (animal) or `P` (phenotypic). `print()` shows the observed trait means, standard
  deviations and pairwise correlations against the replicate distributions, with
  a posterior predictive `p`. Extremes can flag discrepancies, but the `p` is
  not generally uniform under a correct model and a middling value does not
  establish a good fit. Summaries use matching observed-entry masks.

## Breeding values

The animal fit integrates `u` out, so the fit reports `h2` but not what any
individual's additive value is. `breeding_values()` post-processes the fit:

```r
ebv <- breeding_values(fit)
ebv$values        # individuals x traits, response scale
                   # deviations from the fixed effects
ebv$se            # conditional SD, at the posterior means of G, E and beta
ebv$reliability   # 1 - Var(u | observed y) / (A_nn G_kk), in [0, 1]
```

- It is the conditional mean `(A ⊗ G) V^-1 (y - X beta)` with
  `V = A ⊗ G + I ⊗ E`, evaluated at the **posterior means** of `G`, `E` and
  `beta`. That is empirical Bayes: uncertainty in the covariance parameters is
  not propagated, and `se` is conditional on them. It is not a posterior of `u`.
- It needs the `stanimal_data` attribute both fitters attach. `lmm_animal_gibbs()`
  and `lm_phen_cov()` fits do not have one, so they are refused rather than
  silently given a different `A`.
- Missing traits are fine: the gaps are replaced by their posterior mean from
  `imputed_data()` for the point prediction. It is exact for a conditional gap
  mean at fixed plug-in parameters, but these imputations use varying posterior
  draws, so their mean is only an approximation to that quantity. SE and
  reliability condition on the observed entries and include gap uncertainty at
  fixed plug-in parameters, not covariance-parameter uncertainty.
- Use it as a check as well as an output: with `A = I` the prediction reduces to
  shrinking each individual's own residual by `G (G + E)^-1`, and a prediction
  that is not correlated with a simulated `u` means the fit, not the prediction,
  is wrong.
- `se = FALSE` skips the conditional variance. With `m` missing entries the
  SE calculation forms and factorizes a dense `m` by `m` matrix.

## Choosing the model

**Univariate vs multivariate.** `lmm_animal(y, X, A)` estimates a single variance
ratio. `parametrization = "h2"` puts a beta prior on h2 directly, for the
univariate models and for the unstructured multivariate one; `"variance"`
(default) puts half-normals on the two standard deviations and derives h2. Use
`"h2"` when you have a real prior belief about the heritability itself, and when
the data are thin enough that the two parameterisations disagree (section 4). The
ordinary factor model has no switch and always puts `beta_prior` on the trait
heritabilities, for the same reason.

**Residual structure.** `"unstructured"` (default) is a full `E`. `"diagonal"`
assumes independent residuals, which pushes environmental covariance into `G` and
inflates `h2` — the docs warn about this, so prefer `"unstructured"` unless the
design really has independent errors. `"shared_factor"` is a different model, not
a relaxation: the same latent factors drive `G` and `E`, and each factor gets a
heritability.

**`lm_phen_cov()`** is the same machinery without genetics: use it for the
phenotypic covariance, for fixed-effect-adjusted trait covariance, or as a fast
sanity check before an animal model. Its `P` is the **residual** covariance, not
the marginal covariance of `Y`.

**`lmm_animal_laplace()` / `lm_phen_cov_laplace()`** fit the same models by a
Laplace approximation (cmdstanr only): a Gaussian at the mode on the unconstrained
scale, transformed to reported quantities, not a posterior sample. The result has
`$summary` (mean, sd, mad) and `$draws`, no MCMC convergence diagnostics, and the draw-based tools
(`fit_draws()`, `fit_loo()`, `breeding_values()`, `posterior_predictive()`,
`imputed_data()`) refuse it. Standard deviations of variances and heritabilities
are the least reliable entries; use it to explore, then confirm with the sampler.
Arguments in `...` are ignored (with a warning). CmdStan warnings and failures
are surfaced without automatic retries; the seed is recorded in the result.
See `examples/laplace.R` for a sampler comparison.

**`lmm_animal_gibbs()`** is a convenience single-chain sampler with
inverse-Wishart priors. It does **not** reproduce `lmm_animal()`; use it for
comparison or teaching, not as a drop-in.

## Prior choices, hardest first

### 1. `n_factors` — the identification decision

`n_factors` decides how many latent directions exist. Too few underfit the
covariance; too many split one direction across columns and make the factor
summaries multimodal.

- Leaving it `NULL` uses `min(K - 1, max(2, floor(sqrt(K))), 20)` and **warns**,
  because that fallback only sees the number of traits.
- `suggest_n_factors(Y, X)` runs parallel analysis on the residual correlation and
  returns an **upper bound**, not a target. On the package's mouse data it keeps
  7 dimensions while only 2 are identifiable: the 7th clears its cutoff by 0.003
  and the first two carry 28% and 11% of the variance.
- Read the scree (`scree_eigenvalues()`, `plot.scree()`) and stop where the
  eigenvalues flatten. That is the number of factors you can interpret.
- The factors a low-rank model can recover are also bounded: a rank-`k` structure
  cannot reproduce every pairwise correlation. Judge `P` against the best
  rank-`k` approximation, not against zero.

### 2. `column_decay` — the label symmetry

`column_decay` multiplies the prior standard deviation of column `m` by
`column_decay^(m - 1)`. It is the only thing ordering the columns.

- Default `0.8` is mild. `1` means the columns are exchangeable a priori and
  their labels swap between chains.
- **Few factors need more decay than many.** Two columns of comparable size
  needed `0.5` to separate (magnitudes 7.8 / 5.0, `Rhat` 1.0) where `0.8` left
  them at 6.0 / 6.7 with `Rhat` 1.5. But `0.5^(m-1)` shrinks a fifth column to a
  sixteenth of the first, which is wrong for a genuinely many-factor model.

### 3. `hs_global_scale`, `hs_slab_scale` — the loading shrinkage

Regularized horseshoe (Piironen & Vehtari 2017). The slab limits loading prior
standard deviations, not realized loadings. For ordinary factor G and P,
communality is `q / (1 + q)`, where `q` is the random sum of squared raw
loadings after shrinkage. Its prior depends on both scales, the local scales,
factor count and decay, not on the slab alone. The shared-factor model does not
row-normalize its loadings and its communality is on the standardized response
scale. See the reproducible prior simulation in `README.md`; report all the
prior settings when data are weak.

### 4. `lkj_prior` and `beta_prior`

`lkj_prior` (default 2) is used for multivariate unstructured `G` (also with
diagonal `E`) and unstructured `E` (also with factor `G`). It is ignored for
factor `G` with diagonal `E`, shared-factor and univariate models. LKJ(1) is
uniform over correlation matrices, not individual correlations; their
marginals concentrate near zero as the number of traits increases.

`beta_prior` (default `c(1, 1)`) is the prior on h2 for
`parametrization = "h2"`, on the trait heritabilities for the ordinary factor
model, and on the factor heritabilities for
`residual_structure = "shared_factor"`. Both priors are validated **even when the
configuration does not use them**, and supplying an unused one warns instead of
erroring. Do not read a warning about an unused prior as a failure.

The ordinary factor model has no `parametrization` switch, because its two scales
are tied to the phenotypic scale and the prior has to go somewhere: it writes them
as `s_G = sqrt(h2) s_P` and `s_E = sqrt(1 - h2) s_P`, so `beta_prior` is on `h2`
and `s_P` is the free scale. The unstructured multivariate model writes the same
two scales when you ask for `parametrization = "h2"`, and only then. Before 0.7.2
both multivariate models put independent half-normals on `s_G` and `s_E`, which,
because the response is standardized and `s_G^2 + s_E^2` is therefore the
phenotypic variance, is a prior on a circle. The density that induces on h2 is
`1 / sqrt(h2 (1 - h2))`: about 14% of prior mass below 0.05, 14% above 0.95, and a
moderate heritability the least likely answer a priori. A fit at a few hundred
individuals can sit at either end of that range without the data objecting, and
two fits of the same data can then disagree on `h2` and on every `corrG` by an
order of magnitude. If a fit reports an `h2` whose interval touches a boundary,
refit with `parametrization = "h2"` before reading the heritabilities, and treat
the difference between the two fits as the share of the answer the prior was
supplying. No sample-size threshold for that has been validated, so it is a check
and not a rule.

### 5. What is actually identified

| quantity | depends on the orientation? |
|---|---|
| `G`, `E`, `P`, `corrG`, `corrP`, `sigma_G`, `sigma_E`, `h2` | no |
| `communality`, `uniqueness` | no |
| `s_G` (every multivariate fit except the shared-factor one) | no — it is each trait's genetic SD on the standardized response scale, and `communality + uniqueness = 1` makes it `sqrt(G_kk) / y_sd[k]` in a factor fit |
| `loadings`, `factor_magnitude` | **yes.** In the ordinary factor model the likelihood depends on the loadings only through `ΛΛ'`, so it is flat along every rotation and only the prior picks a point |
| `factor_heritability` | **partly.** The shared-factor likelihood is invariant to sign flips and column permutations, not to general rotations |

Two solutions with the same `G` can differ by a rotation, a sign flip or a
permutation of the columns. So a high `Rhat` on `loadings`/`factor_magnitude` with
clean `G`/`communality` is an orientation that moved, not a broken fit. A high
`Rhat` on `communality` **is** a real problem: it means genuinely different
decompositions were visited, so reduce `n_factors`.

**Sign flips are handled for you.** `fit_draws()` — and therefore `draws_mean()`
and `draws_quantile()` — flips each draw so that the trait with the largest mean
absolute loading on a factor is positive, chosen once over the whole posterior.
Without that, a chain that visited both orientations averages to zero and the
factor silently disappears from a mean or a heatmap. `fit$draws("loadings")` is
still the raw, unaligned output. Signs are the only thing repaired: a permutation
of the columns, or a genuine rotation between draws, is not, which is why those
are judged on `G` and `communality`.

Triage when a factor fit looks bad, in order: check `P`/`G` and `communality`
`Rhat` → if those are clean, ignore the loadings → if not, lower `n_factors` →
then lower `column_decay`.

### 6. Checking a prior

```r
plot_prior_posterior(fit)   # compares with the prior the fit recorded
animal_prior_draws(3)       # standalone prior simulator for the unstructured model
```

Both `plot_prior_posterior()` and `plot_forest()` reconstruct the prior for
unstructured and ordinary factor-G animal fits with complete recorded metadata.
For factor G, they simulate the fitted horseshoe and column decay; `lkj_prior`
does not set genetic correlations. Shared-factor or custom fits, bare draws
and older fits missing metadata need explicit `prior` draws in the
format returned by `animal_prior_draws()`, with matching trait names and
parameter types; the caller must ensure the distribution is the model's prior.
Supplying a different `lkj_prior` warns. For example, for bare draws from a known
unstructured variance model, use `trait_names = c("a", "b")` and
`prior = animal_prior_draws(2, trait_names = c("a", "b"))`. For factor models,
simulation with `random_covariance(method = "factor")` builds a known
covariance, and recovery of `G`/`P` can be scored directly.

## Anti-patterns

- **Do not fit a factor model on a rotated basis.** After `size_shape()` every
  loading is a dense mixture, so the sparsity prior has nothing to shrink and the
  low-rank structure sits in the rotation instead of in the model.
- **Do not read `factor_magnitude` as "importance"** unless the orientation is
  pinned by `column_decay`.
- **Do not compare covariances with an absolute RMSE.** Entries can be ~1e-3, so a
  small absolute number can be a large relative error. Use the correlation scale
  and a relative error.
- **Do not assume `P` is the marginal covariance.** It is the covariance of the
  residuals after `X`.
- **Do not assume chains share a start.** Each chain draws its own start,
  uniform on (-1, 1) on the unconstrained scale (`init = 1`, narrower than the
  engines' (-2, 2), which stalls correlation factors over many traits). Pass
  `init` if you need reproducibility at a specific point.

## Where to look next

- `README.md` — the full narrative, worked examples, and the mathematics of each
  model.
- `examples/` — runnable scripts; `examples/README.md` says what each covers.
- `?lmm_animal`, `?lm_phen_cov`, `?size_shape`, `?makeA` — argument-level detail.
- `DEVELOPMENT.md` — **only if you are changing stanimal itself**: repository
  layout, build and test commands, and the packaging conventions.
