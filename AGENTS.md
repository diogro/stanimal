# Guide for agents working with stanimal

`stanimal` fits Bayesian quantitative-genetic models in Stan (CmdStan via cmdstanr,
with an rstan fallback). Everything is returned on the **original data scale**, so
estimates can be compared with the trait values directly.

This file is the interface guide: what to call, what comes back, and how to make
the prior choices that the function signatures do not make for you. Usage
examples live in `README.md` and `examples/`; the package sources are in
`package/`.

## Quick decision table

| Goal | Call |
|---|---|
| One trait, heritability | `lmm_animal(y, X, A)` |
| Several traits, covariance structure | `lmm_animal(Y, X, A)` |
| Many traits, `G` needs to be sparse/low-rank | `lmm_animal(Y, X, A, genetic_structure = "factor")` |
| `G` and `E` should share the same latent factors | add `residual_structure = "shared_factor"` (needs `"factor"`) |
| No pedigree, only the phenotypic covariance | `lm_phen_cov(Y, X)` |
| Phenotypic covariance of many traits | `lm_phen_cov(Y, X, covariance_structure = "factor")` |
| Build `Y`/`X` from a formula, with `A` | `animal_design(formula, data, pedigree = ...)` |
| Pedigree with related founders/metafounders | `makeA(pedigree, Gamma = , groups = )` |
| Markers instead of a pedigree | `genomic_relationship(M)` → pass as `A` |
| How many factors? | `suggest_n_factors(Y, X)` + `scree_eigenvalues(Y, X)` |
| Many correlated traits on one size gradient | `size_shape()` then `back_transform()` — but only around an **unstructured** model, see the anti-patterns |

`A` is the third positional argument of `lmm_animal()`, so `lmm_animal(Y, X, A)`.

## The interface

```r
lmm_animal(Y, X = NULL, A, lkj_prior = 2, beta_prior = c(1, 1),
           parametrization = c("variance", "h2"),
           iter = 3000, warmup = floor(iter/2), chains = 4, cores = chains,
           model_file = NULL, control = list(adapt_delta = 0.99),
           cmd_stan = TRUE, blend_A = FALSE,
           residual_structure = c("unstructured", "diagonal", "shared_factor"),
           genetic_structure = c("unstructured", "factor"), n_factors = NULL,
           hs_global_scale = 1, hs_slab_scale = 2, column_decay = 0.8, ...)

lm_phen_cov(Y, X = NULL, lkj_prior = 2, iter = 3000, warmup = floor(iter/2),
            chains = 4, cores = chains, control = list(adapt_delta = 0.99),
            cmd_stan = TRUE, model_file = NULL,
            covariance_structure = c("unstructured", "factor"), n_factors = NULL,
            hs_global_scale = 1, hs_slab_scale = 2, column_decay = 0.8, ...)

lmm_animal_gibbs(Y, X = NULL, A, iter = 1000, warmup = floor(iter/2), thin = 1)
```

`Y` is a numeric matrix (individuals × traits) or a vector; `X` defaults to an
intercept; `A` is a symmetric positive-definite relationship matrix whose rows
are matched to `Y` by name when both have names. `...` is forwarded to the
sampler, so `seed`, `refresh`, `init` and `show_messages` work there.

Both model functions return a `CmdStanMCMC` (or a `stanfit` with
`cmd_stan = FALSE`) carrying the configuration they ran with:
`attr(fit, "genetic_structure")`, `"residual_structure"`, `"parametrization"`,
`"lkj_prior"`, `"beta_prior"`, `"n_factors"`, `"x_names"`, `"y_names"`.

### What each configuration reports

| fit | quantities |
|---|---|
| univariate | `sigma_G`, `sigma_E`, `h2`, `beta` |
| multivariate, unstructured | `G`, `E`, `P`, `corrG`, `sigma_G`, `sigma_E`, `h2`, `beta`, `log_lik` |
| factor `G` | + `loadings`, `uniqueness`, `communality`, `factor_magnitude`, `s_G` |
| shared factor | + `loadings`, `uniqueness`, `communality`, `factor_magnitude`, `factor_heritability`, `genetic_uniqueness`, `residual_uniqueness` |
| `lm_phen_cov()` | `P`, `corrP`, `sigma_P`, `beta`, `log_lik` (+ `loadings`, `uniqueness`, `communality`, `factor_magnitude` for `"factor"`) |

`h2` exists only for animal fits. `log_lik` is per individual, for LOO/WAIC, and
only the multivariate and `lm_phen_cov()` models report it — univariate fits do
not, so LOO cannot be computed for them.

### Reading the draws

```r
fit$summary(c("G", "h2"))          # raw variables, tidy data frame
fit_draws(fit)                     # list of arrays: traits first, draws last
draws_mean(fit, "G")               # posterior mean as a named matrix
draws_quantile(fit, "corrG")       # 2.5% and 97.5% as a 3-d array
rename_beta(fit)                   # beta named after the columns of X
to_draws_list(fit)                 # one posterior draws matrix per quantity
```

`draws_mean()`/`draws_quantile()` take `parameter = NULL` for everything, and a
`traits` argument to relabel. Factor extras keep their natural shape:
`loadings` is **trait × factor × draw**, and `factor_magnitude` and
`factor_heritability` are **factor × draw**.

## Choosing the model

**Univariate vs multivariate.** `lmm_animal(y, X, A)` estimates a single variance
ratio. `parametrization = "h2"` puts a beta prior on h2 directly (univariate
only); `"variance"` (default) puts half-normals on the two standard deviations and
derives h2. Use `"h2"` when you have a real prior belief about the heritability
itself.

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

Regularized horseshoe (Piironen & Vehtari 2017). `hs_global_scale` (default 1) is
the half-Cauchy scale on the global shrinkage and sets roughly how many loadings
survive; `hs_slab_scale` (default 2) bounds the largest loading. Larger
`hs_global_scale` → denser factors. For a scale derived from an expected number
of nonzero loadings, see the reference in `?lmm_animal`.

### 4. `lkj_prior` and `beta_prior`

`lkj_prior` (default 2) shrinks correlations toward zero and is used only by the
multivariate unstructured `E` (including factor `G` with unstructured `E`). It is
**ignored** for `"diagonal"` and `"shared_factor"`, and explicitly passing it
there warns.

`beta_prior` (default `c(1, 1)`) is the prior on h2 for
`parametrization = "h2"` and on the factor heritabilities for
`residual_structure = "shared_factor"`. Both priors are validated **even when the
configuration does not use them**, and supplying an unused one warns instead of
erroring. Do not read a warning about an unused prior as a failure.

### 5. What is actually identified

| quantity | depends on the orientation? |
|---|---|
| `G`, `E`, `P`, `corrG`, `corrP`, `sigma_G`, `sigma_E`, `h2` | no |
| `communality`, `uniqueness` | no |
| `s_G` (ordinary factor model) | no — it is each trait's genetic SD on the standardized response scale, and `communality + uniqueness = 1` makes it `sqrt(G_kk) / y_sd[k]` |
| `loadings`, `factor_magnitude` | **yes.** In the ordinary factor model the likelihood depends on the loadings only through `ΛΛ'`, so it is flat along every rotation and only the prior picks a point |
| `factor_heritability` | **partly.** The shared-factor likelihood is invariant to sign flips and column permutations, not to general rotations |

Two solutions with the same `G` can differ by a rotation, a sign flip or a
permutation of the columns. So a high `Rhat` on `loadings`/`factor_magnitude` with
clean `G`/`communality` is an orientation that moved, not a broken fit. A high
`Rhat` on `communality` **is** a real problem: it means genuinely different
decompositions were visited, so reduce `n_factors`.

Triage when a factor fit looks bad, in order: check `P`/`G` and `communality`
`Rhat` → if those are clean, ignore the loadings → if not, lower `n_factors` →
then lower `column_decay`.

### 6. Checking a prior

```r
plot_prior_posterior(fit)   # compares with the prior the fit recorded
animal_prior_draws(3)       # standalone prior simulator for the unstructured model
```

A sparse factor fit has no closed-form prior to redraw, so `plot_prior_posterior()`
errors there by design; pass `prior` explicitly if you need that plot. The
reliable check for the factor models is simulation: `random_covariance(method =
"factor")` builds a known covariance, and recovery of `G`/`P` can be scored
directly.

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
- **Do not assume chains share a start.** No `init` is forced: each chain draws
  its own. Pass `init` if you need reproducibility at a specific point.

## Working in this repository

- `package/` is the R package; `README.md` and `examples/` are outside it.
- Test: `Rscript -e 'testthat::test_local("package", filter="model-P", reporter="summary")'`
- Docs: `cd package && Rscript -e 'roxygen2::roxygenise(".")'`
- Check: `R CMD build package && R CMD check --no-tests --no-manual --no-build-vignettes stanimal_0.6.0.tar.gz`
  (checking the directory directly fails: DESCRIPTION has `Authors@R` but no
  generated `Author:` field)
- Examples call `library(stanimal)`, so reinstall (`R CMD INSTALL package`) after
  changing anything under `package/R` or `package/inst/stan`.
- `test-readme.R` executes the root README's R blocks but **skips inside the built
  tarball**, so a passing `R CMD check` says nothing about the README. Run it from
  a checkout. Its blocks share one environment, so a snippet must not clobber a
  name a later one uses.
- Roxygen markdown is **off**: use `\code{}`/`\strong{}` and escape prose `%` as
  `\%`, or the rest of the Rd line is silently commented out.
- `.Rbuildignore` controls the tarball; `.gitignore` does not. Compiled Stan
  executables must stay in `.Rbuildignore`.
- Fit tolerances are calibrated: measure the actual error, then set roughly twice
  it. Loose absolute bounds have hidden real bugs here before.
