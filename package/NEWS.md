# stanimal (development version)

* The unstructured `lm_phen_cov()` model no longer multiplies out its Cholesky
  factor and factorizes the product again for each pattern of observed traits.
  A pattern that observes every trait uses the factor directly; an incomplete
  one takes it from a QR decomposition of the factor's rows. Refactorizing
  failed in floating point when the factor had tiny diagonal entries, so a
  35-trait fit could reject every initial value. This regressed in 0.8.0 with
  incomplete trait records, and complete-data fits also save one `K x K`
  Cholesky per gradient.

* The default `adapt_delta` drops from 0.99 to 0.9 for `lm_phen_cov()` and for
  every `lmm_animal()` structure except `residual_structure = "shared_factor"`,
  which keeps 0.99. The lower target takes larger steps: in the examples,
  factor and structure-comparison fits run 2 to 3 times faster with no
  divergences. The shared-factor fit diverged at 0.9 in 6 of 10 seeds, so it
  does not. `lmm_animal()`'s `control` now defaults to `NULL`, which picks the
  target from the structure. Pass `control = list(adapt_delta = 0.99)` to get
  the old behaviour.

* Laplace fits (`lmm_animal_laplace()`, `lm_phen_cov_laplace()`) label their
  summary and draws with the column names of `Y` and `X` (`G[wing,tail]`,
  `h2[tail]`, `beta[tail,size]`), as `fit_draws()` names sampler fits, and
  record them as `y_names` and `x_names`. An unnamed `Y` keeps Stan's indices.
  Code that looked quantities up as `G[1,1]` on a named `Y` must use the
  names.

# stanimal 0.8.1

## Bug fixes

* Two fits of the same model in one session that started from the same random
  state within a minute (for example with `set.seed()` before each fit) wrote
  to the same CmdStan output files. The first fit then silently returned the
  second fit's draws, and garbage-collecting the first fit deleted the
  second's files. cmdstanr builds output names from the minute, the process id
  and one `runif()`; stanimal now gives every sampler run a unique
  `output_basename` (from `tempfile()`, which `set.seed()` does not affect).
  Laplace fits read their draws at once and keep no fit, so they were not
  affected.
  An `output_dir` or `output_basename` passed through `...` is respected.

# stanimal 0.8.0

## New features

* `lmm_animal()` gains `eigen_A`, a precomputed eigendecomposition of `A`
  (`list(Q =, lambda =)`, e.g. from `eigen(A, symmetric = TRUE)`). A caller
  fitting many models on the same relationship matrix, such as one model per
  trait subset, pays the `O(N^3)` decomposition once instead of once per fit.
  It works with every structure and parametrization, is checked against `A`
  (reconstruction and orthonormality) before fitting, follows `A` when it is
  reordered to match `Y`, and takes `blend_A`. The decomposition can come from
  any solver, including a GPU one.

## Changes in behaviour

* The animal models take the eigendecomposition of `A` as data instead of
  decomposing `A` themselves in Stan. `lmm_animal()` and
  `lmm_animal_laplace()` decompose it once per fit in R; previously every
  chain repeated it, computing the eigenvectors and the eigenvalues
  separately. The model is unchanged: the eigenvalues are passed in Stan's
  ascending order, and in the tests a model that decomposes `A` the old way
  reproduces the new one draw for draw. Where `A` has repeated eigenvalues
  the two libraries may choose different bases for them, which can change
  the MCMC trajectory through rounding but not the posterior. A custom
  `model_file` that declares `A` still receives it.
* The decomposition now runs on R's BLAS/LAPACK. With an optimised library
  (OpenBLAS, MKL, Accelerate) it is several times faster than Stan's; with
  R's reference LAPACK it is about 1.5 times slower, so a single large fit
  can take slightly longer in wall clock (about 6% at N = 4000) while using a
  third of the CPU. Use an optimised BLAS for large `N`, or pass `eigen_A`.
* `fit_loo()` and the `loo()` fit methods now refuse LOO and WAIC for animal
  fits with missing responses; fits without response metadata warn that their
  likelihood terms cannot be verified.
* `plot_G_ellipses()` now labels each trait pair with its own genetic
  correlation and interval when there are four or more traits.

## Bug fixes

* Animal posterior predictive replicates now use the correct relationship
  covariance between individuals.
* `breeding_values()` retains imputation-based predictions but now computes
  standard errors and reliabilities conditional only on observed responses.
* Univariate breeding values use the posterior mean variances, not squared
  posterior mean standard deviations.
* Predictive summaries now use the same observed-entry masks as the data.
* `size_shape()` now reports an actionable error when the requested shape axes
  exceed the numerical rank of the centered shape data.
* `genomic_relationship()` excludes fully missing markers even with supplied
  frequencies and `max_missing = 1`, warning when they would otherwise pass
  the filters.
* Relabeling `fit_draws` now also updates the trait labels of factor loadings.
* Missing-data EM starts with a positive-definite covariance and errors on
  non-convergence or singularity; design diagnostics warn if unavailable.

## Documentation

* Corrected the opossum data help and example and agent guides.
* Corrected LKJ-use warnings, prior and model-comparison guidance, and
  posterior predictive tail-probability explanations.

# stanimal 0.7.4

## Bug fixes

* `breeding_values()` now computes reliability using the individual's prior
  variance, `A[n,n] * G[k,k]`, instead of assuming unit diagonals in `A`.

* Laplace fits write the optimization mode at high precision, avoiding a
  Cholesky round-trip failure caused by insufficient CSV precision. They no
  longer suppress warnings or retry unrelated errors; a failed Laplace stage
  reports CmdStan output. Explicitly supplied unused model arguments warn as in
  the samplers.

## Changes in behaviour

* `plot_prior_posterior()` and `plot_forest()` can now reconstruct the ordinary
  factor-G prior by simulating its recorded horseshoe and factor settings.
  Shared-factor fits still require explicitly supplied prior draws.

* `plot_prior_posterior()` and `plot_forest()` now require explicit prior draws
  for bare draws or fits missing the metadata needed to reconstruct the prior,
  rather than silently assuming the default model. They validate supplied prior
  labels and use recorded trait names, and warn when an LKJ override differs
  from the fitted prior.

## Documentation and tests

* Added a runnable Laplace example and tests covering model configurations,
  missing responses, error reporting, and a comparison with sampling.
* Clarified prior, breeding-value and data documentation. Removed committed
  build and check artifacts and updated the development workflow to keep them
  outside the checkout.

# stanimal 0.7.3

## Documentation

* `ratones_morpho` and `ratones_pedigree` stay, with a weak identification
  section explaining what they are: a long-term artificial-selection experiment
  that happens to carry a pedigree, in which 329 individuals from 58 sires
  identify a 35-trait `G` only weakly, so a heritability read from
  `lmm_animal()` on these traits depends on the prior as much as on the data.
  They remain the small worked example; `opossum` and `opossum_grm` are the
  datasets for anything about genetic covariance.

## New data

* `opossum` and `opossum_grm` ship with the package: 553 Monodelphis domestica
  from Porto et al. (2016), measured on the same 35 cranial landmarks as
  `ratones_morpho`, together with a genomic relationship matrix built from their
  3,696 markers. The data are stored untransformed, as the published analysis
  used them. The genomic matrix is the right relationship matrix for these
  animals: the colony pedigree shipped with the archive puts every pair at a
  kinship of 0.69 to 0.93, which leaves additive variance indistinguishable from
  a shared colony effect.

## Changes in behaviour

* `hs_slab_scale` now defaults to 1 instead of 2. The earlier explanation that
  the slab alone determines expected communality was incorrect: the entire
  horseshoe prior, factor count and decay determine its distribution. With 600 individuals and 12 traits
  the factor model recovers a simulated `G` to a relative MSE of 0.07 at a slab
  of 2, while on the mice data the slab moved the posterior mean heritability
  through 0.10, 0.21, 0.78 and 0.80 at slabs of 0.5, 1, 2 and 4. A factor fit is
  therefore reproducible only if the slab is reported with it, and
  `hs_slab_scale` should be reported alongside the other factor prior settings.
  `hs_global_scale` is left at 1; across 0.1 to 4 the same mouse fit moved its
  heritability between 0.785 and 0.790, a result specific to that example.

# stanimal 0.7.2

## New features

* `parametrization = "h2"` now applies to the unstructured multivariate model as
  well as to the univariate ones. It writes the two scales as a phenotypic scale
  and a heritability per trait, `s_G = sqrt(h2) s_P` and `s_E = sqrt(1-h2) s_P`,
  and puts `beta_prior` on the trait heritabilities. The default is unchanged
  (`"variance"`), and the two parameterizations agree when the data are
  informative, so this is the switch to reach for when there are few individuals
  per trait or when the heritability itself is the object of inference.

## Bug fixes

* The ordinary factor model now puts the beta prior on the trait heritabilities
  instead of on two standard deviations. The response is standardized trait by
  trait, so `s_G^2 + s_E^2` is the phenotypic variance and the data pin it near
  1, which made the old half-normal priors on `s_G` and `s_E` a prior on a
  circle: the density induced on h2 was `Beta(1/2, 1/2)`, with about 14% of its
  mass below 0.05 and 14% above 0.95 and a moderate heritability the least likely
  answer a priori. At a few hundred individuals the likelihood is weak enough for
  that to show through, and a fit could sit at either end of the range, reporting
  h2 near 0 or near 1 with genetic correlations to match. The model is now
  parameterized by a phenotypic scale `s_P` and a heritability `h2`, with
  `s_G = sqrt(h2) s_P` and `s_E = sqrt(1 - h2) s_P`, so `beta_prior` (default
  `c(1, 1)`, flat on h2) is the prior on the trait heritabilities and `s_P` is
  the free scale. `h2` is still reported and still equals `diag(G) / diag(P)`.
  The fits record `parametrization = "h2"`, and `s_P` is a reported parameter.

# stanimal 0.7.1

## New features

* `lmm_animal_laplace()` and `lm_phen_cov_laplace()` fit the same models by a
  Laplace approximation at the posterior mode instead of sampling them. They
  return posterior means and standard deviations for the covariance quantities
  and finish in a fraction of the time, and they deliberately return no draws:
  there are no chains, so there are no convergence diagnostics, and the
  draw-based functions refuse the result. The standard deviation of a variance
  or a heritability is the least reliable number they report, because a Gaussian
  at the mode misstates a skewed posterior.

# stanimal 0.7.0

## New features

* Incomplete trait records are used instead of discarded. `Y` may contain `NA`
  in `lmm_animal()` and `lm_phen_cov()`: a record is kept as long as it has one
  observed trait, and the missing values no longer require dropping the row and
  losing the traits that were measured. `lm_phen_cov()` evaluates each
  missingness pattern on its observed subvector of `P`, so the gaps are
  integrated out and cost no parameters; `lmm_animal()` fills them with
  parameters (`Y_mis`, reported in the draws), which is the exact marginal over
  the observed data and lets relatives inform a gap. A fully missing record is
  allowed in the animal model and an error elsewhere.
* `animal_design()` no longer drops rows that miss only the response; it still
  drops rows with missing predictors, and `na.action` now applies to those.
* `size_shape()` accepts a response with gaps, which it could not before: size is
  a mean of the traits and shape is a difference between them, so neither is
  defined for a row with a missing value. The response is completed first, by a
  maximum-likelihood multivariate normal estimated from the values that were
  measured (EM), and each gap becomes its conditional mean given the traits that
  individual has. The rest of the transformation is unchanged. `$gaps` reports
  how many entries were completed and `$scores_sd` the standard deviation of each
  completed score, which is 0 for the individuals that were measured in full.
  This is a point completion, so a fit on the resulting scores treats them as
  data: the spread in `$scores_sd` does not travel into that posterior. A row
  with nothing observed, and a pair of traits never measured together, are
  errors. `pc1_variance_share()`, `plot_G_ellipses()` and the diagnostic
  `animal_design()` runs on its response accept gaps for the same reason.
* `imputed_data()` returns completed datasets, one per posterior draw: the
  observed entries as they are and the gaps filled on the response scale. An
  animal fit's gaps are the parameters it sampled; a phenotypic fit has none, so
  its gaps are drawn from the normal conditional distribution given that row's
  observed traits. `print()` reports each gap with its 95% interval. These are
  multiple imputations, not a single fill, and `fit$draws("Y_mis")` alone is not
  enough to reconstruct them: it is on the standardized scale and says nothing
  about which individual and trait each value belongs to.
* `lm_phen_cov()` reports pointwise `log_lik` for the unstructured and factor
  fits, so `loo` and `waic` work for every model in the package. It is the
  density on the standardized response, as in the animal models.
* Univariate fits report pointwise `log_lik`, so every model in the package now
  exposes the per-record terms that `loo` and `waic` need. The h2 and variance
  parametrizations describe the same model and report the same density.
* `fit_loo()` compares fits by leave-one-out cross-validation, wrapping the `loo`
  package on the `log_lik` draws with the chain structure kept; `criterion =
  "waic"` gives WAIC, and `loo(fit)` dispatches to it once `loo` is loaded. The
  terms are the model block's, one per record of the rotated response, so this
  compares model structures on the same data rather than leaving out one
  individual. `loo` is a new Suggests.
* `posterior_predictive()` replicates the response from the posterior using the
  reported `G`, `E` or `P`, and reports the observed trait means, standard
  deviations and pairwise correlations against the replicates, with a tail
  probability for each (`print()`) and the distributions themselves (`plot()`).
  This is the individual-level check that `fit_loo()` is not.
* `breeding_values()` predicts each individual's additive genetic value, which
  the animal models integrate out and never report. It is the conditional mean
  `(A %x% G) V^-1 (y - X beta)` at the posterior means of `G`, `E` and `beta`,
  computed per individual through the eigendecomposition of `A`, and returns the
  predictions on the response scale with an optional standard error and the
  reliability `1 - Var(u | y) / diag(G)`. Being empirical Bayes, it does not
  propagate the uncertainty in the covariance parameters. Missing values are
  filled by their conditional mean first, which is exact for this linear
  predictor. Fits from `lmm_animal()` carry the data it needs; `lm_phen_cov()`
  and `lmm_animal_gibbs()` do not.

## Changes

* `fit_draws()`, `draws_mean()` and `draws_quantile()` now sign-align the factor
  loadings. A loading column is identified only up to its sign, so a chain can
  visit both orientations, and the average of the raw draws then cancels the
  factor out of a posterior mean or a heatmap. Each draw is flipped so that the
  trait with the largest mean absolute loading is positive, chosen once over the
  whole posterior. `fit$draws("loadings")` still returns the draws as the sampler
  wrote them; nothing else is reoriented, since column permutations and rotations
  are a different problem and are judged on `G` and `communality`.

## Documentation

* `AGENTS.md` states what each configuration reports, how the priors behave, and
  which quantities a factor model can and cannot identify.

# stanimal 0.6.0

## New features

* `lmm_animal()` gains `genetic_structure = "factor"`, which models the genetic
  covariance as a sparse factor structure on the correlation scale,
  `G = diag(s_G) (Lambda Lambda' + Psi) diag(s_G)`, with regularized-horseshoe
  priors on the loadings (Piironen and Vehtari 2017). `n_factors` sets the
  number of factors; if it is omitted, the factor model starts from the
  sublinear default `min(K - 1, max(2, floor(sqrt(K))), 20)` and reports the
  value. `hs_global_scale` and `hs_slab_scale` control the shrinkage, and
  `column_decay` optionally multiplies the prior standard deviation of factor
  column `m` by `column_decay^(m - 1)` to favour decreasing factor magnitudes.
  The unstructured model remains the default. Factor fits report
  `loadings`, `factor_magnitude`, `uniqueness`, `communality` and `s_G` in
  addition to `G`, `E`, `P`, `h2` and the other summaries.
* `lmm_animal()` gains `residual_structure = "shared_factor"` for the
  shared-factor model
  `G = Lambda diag(h2_factor) Lambda' + diag(psi_G)` and
  `E = Lambda diag(1 - h2_factor) Lambda' + diag(psi_E)`.
  The same latent factors structure both covariances, and the factor
  heritabilities identify how much of each factor is genetic; their beta prior
  is controlled by the existing `beta_prior` argument. It requires
  `genetic_structure = "factor"` and reports `loadings`, `factor_magnitude`,
  `factor_heritability`, `genetic_uniqueness`, `residual_uniqueness`,
  `uniqueness` and `communality`.
* `lm_phen_cov()` gains `covariance_structure = "factor"`, which replaces the
  unstructured phenotypic covariance with a sparse factor model on the
  correlation scale. It uses the same regularized-horseshoe loadings,
  `n_factors` default, `hs_global_scale`, `hs_slab_scale` and `column_decay` as
  the animal factor model, and reports `loadings`, `factor_magnitude`,
  `uniqueness` and `communality`. The factor fit estimates the residual
  standard deviations, so `P`, `corrP` and `sigma_P` describe the covariance
  left after the fixed effects rather than the marginal covariance of `Y`.
* `suggest_n_factors()` chooses the number of factors from the data by parallel
  analysis (Horn 1965), comparing the eigenvalues of the residual trait
  correlation matrix against those of data with no shared structure. The result
  is an upper bound, and its attributes carry the full eigenvalue comparison.
  `scree_eigenvalues()` returns the eigenvalues behind a scree plot, and
  `plot()` draws it with the parallel-analysis cutoffs. Both residualize on the
  fixed effects first, so they describe the correlation the factor model fits.
* The multivariate models now report pointwise `log_lik`, enabling LOO, WAIC
  and other predictive comparisons with the fits.
* `fit_draws()` handles the factor outputs: factor fits return `loadings` as a
  trait by factor by draw array, `uniqueness` and `communality` as trait by
  draw matrices, and `factor_heritability` and `factor_magnitude` as factor by
  draw matrices.
* Prior arguments are validated consistently even when the selected model does
  not use them, and explicitly supplying an unused prior now warns instead of
  being silently ignored.
* `examples/factor_model.R` demonstrates the sparse factor model on a
  twelve-trait modular hierarchy, and `examples/structure_comparison.R`
  compares the unstructured, diagonal, factor and shared-factor covariance
  structures.

## Bug fixes

* Fits from the `rstan` backend now keep the trait and term names, so
  `fit_draws()` labels its arrays with the columns of `Y` and `rename_beta()`
  with the columns of `X`. Those attributes were only attached on the cmdstanr
  path, which left `rstan` fits with generic `trait1`, `trait2` labels.
* `plot_prior_posterior()` and `plot_forest()` compare against the prior the fit
  actually used instead of always drawing the unstructured one. `lmm_animal()`
  now records `lkj_prior`, `beta_prior` and `parametrization` on the fit, so an
  unstructured fit is compared with its own LKJ shape and a univariate
  `parametrization = "h2"` fit with its beta prior on the heritability (which
  `animal_prior_draws()` can now generate through `beta_prior`). A sparse factor
  fit has no closed-form prior to reproduce, so that combination errors and asks
  for an explicit `prior` rather than quietly plotting the wrong comparison.
* The factor phenotypic model in `lm_phen_cov()` no longer pins the diagonal of
  `P` to the marginal variance of `Y`. Removing the fixed effects lowers the
  residual variance, so the model could only fit the likelihood by inflating
  the correlations. It now estimates a residual scale per trait, as the
  unstructured model does, and `sigma_P` is the residual standard deviation.

## Changes

* `column_decay` now defaults to `0.8` instead of `1` in both factor models. With
  no decay the factor columns are exchangeable a priori, so their labels are free
  to swap between chains and `factor_magnitude` reports large `Rhat` values on a
  healthy fit; the mild decay orders the columns without forcing a hierarchy.
  `column_decay = 1` restores the previous behaviour. Supplying `column_decay`
  without a factor model now warns instead of erroring, like the other factor
  priors.
* Leaving `n_factors` to the fallback now warns with the number it used, instead
  of reporting it with a message. The fallback only sees the number of traits;
  `suggest_n_factors()` is the data-driven answer.
* Chains no longer all start at the same point. The sampler used to be called
  with `init = 0`; each chain now draws its own start from the engine default.
  A shared start makes chains agree by construction, so `Rhat` cannot detect a
  posterior mode that no chain visited. Pass `init` through `...` to pin the
  start of a run that is hard to initialize.

# stanimal 0.5.1

## Name change

* The package is now called **stanimal** instead of **stanAnimal**. The
  functions keep their names (`lmm_animal()`, `animal_design()`, `makeA()`), so
  only `library()`, the installed package name and the two cache settings
  change: `options(stanAnimal.cache_dir = )` is now
  `options(stanimal.cache_dir = )`, and `STANANIMAL_CACHE_DIR` is now
  `STANIMAL_CACHE_DIR`. Compiled models are cached per package name, so the
  first fit after upgrading recompiles; `stanimal_clear_cache()` removes the old
  files.

## New features

* `makeA()` and `animal_design()` accept `Gamma`, a relationship matrix among
  founders, and `groups`, an optional founder-to-group assignment. This builds
  the `A^Gamma` relationship matrix for related founders and metafounders
  (Legarra et al. 2015): the founder block is initialized from `Gamma` instead
  of the identity, and the tabular recursion propagates those relationships.
  The default remains the standard pedigree with unrelated founders.
* `random_covariance()` draws random covariance matrices for simulations, either
  from the LKJ distribution used as a prior by `lmm_animal()` (default), from
  the partial-correlation construction of `evolqg::RandomMatrix()`, or from a
  positive size factor mixed with a sign-oriented random LKJ residual
  (`method = "factor"`). The factor method solves the factor strength (or, for
  weak integration, a shrinkage toward independence) to match a requested
  integration level, keeps PC1 positive, and redraws the residual for every
  matrix, so draws at a fixed integration still differ. It spans the full range
  without the near-singular conditioning of bending an LKJ spectrum directly.
  The port handles one and two traits, where the original indexing fails, and
  both the LKJ and factor methods avoid the nearly singular matrices the
  partial method often produces. `trait_comparison.R` uses it for the
  simulated truth.
* `fit_draws()` extracts the draws of a fit as arrays indexed by trait, with
  draws last and every dimension named: `G`, `E`, `P`, `corrG`, `corrE`,
  `corrP` are trait by trait by draw arrays and `h2`, `sigma_G`, `sigma_E` and
  `sigma_P` are trait by draw matrices, so posterior means and correlations no
  longer need a helper that matches `G[i,j]` variable names. It accepts the
  fits of `lmm_animal()`, `lm_phen_cov()` and `lmm_animal_gibbs()` (whose `R` is
  returned as `E`), `posterior` draws objects and `back_transform()` output,
  and derives quantities the fit does not report (`P`, `corrE`, `corrP`,
  `sigma_G`, `sigma_E`, `sigma_P`, `h2`).
* `draws_mean()` and `draws_quantile()` summarise a fit or an `fit_draws()`
  object without the draw dimension, so `draws_mean(fit, "G")` is the posterior
  mean as a trait by trait matrix and `draws_quantile(fit, "corrG")` gives its
  interval. The examples and the README use them instead of defining their own
  `post_mean()` helper.
* `lmm_animal()` and `lm_phen_cov()` keep the response column names on the
  cmdstanr fit, so `fit_draws()` and the summary functions label the trait
  dimensions without an extra argument.

## Improvements

* Declare the Stan project's R-universe as an additional repository for
  `cmdstanr` and use it in the installation guidance.
* Update package metadata and DOI references for the CRAN release.

* `lmm_animal_gibbs()` validates its arguments with the same helpers as the Stan
  models instead of its own inline checks: `Y` and `X` through the shared
  response and design check, `A` through the shared symmetry and definiteness
  checks, and `iter`, `warmup` and `thin` through the shared iteration check.
  Its messages now name the argument the caller passed.

* Input checks are uniform across the package: every error names its argument in
  backticks and is raised without a call, so the message points at `Y`,
  `lkj_prior` or `maf` rather than at the internal function that noticed. The
  checks live in `utils.R` and are reused, instead of each function testing its
  own conditions inline.

* `lmm_animal()` and `lm_phen_cov()` check `iter`, `warmup`, `chains`, `cores`,
  `cmd_stan` and, for `lmm_animal()`, `blend_A` before fitting, so a warmup
  longer than the run is an error naming the argument instead of a failure
  inside the Stan call.

* `X` is optional in `lmm_animal()`, `lm_phen_cov()` and `lmm_animal_gibbs()`.
  Leaving it out fits an intercept-only model: one column of ones named
  `(Intercept)`, so the draws come out as `beta[(Intercept)]`, as they do with
  an explicit intercept column, and `beta` keeps its trait dimension.

# stanAnimal 0.5.0

## Breaking changes

* Three functions are renamed, with no deprecated aliases, so a script that
  uses the old names fails at the call:
  * `animal_draws()` is now `rename_beta()`, which says what it does: it takes
    the draws from a fit and renames the `beta` variables after the columns of
    `X`. The old name read as the way to get draws at all, which is
    `fit$draws()`.
  * `pc_dominance()` is now `pc1_variance_share()`, which names the first-PC
    variance proportion without suggesting a genetic dominance effect.
  * `lmm_multiCov()` is now `lm_phen_cov()`: the phenotypic covariance model
    has no random effects, so it is a linear model rather than a mixed model.

* `back_transform()` returns a list of class `size_shape_draws` instead of a
  single `draws_matrix`: `G`, `E`, `P`, `corrG` and `corrP` are trait by trait
  by draw arrays and `h2` is a trait by draw matrix, the same layout as
  `lmm_animal_gibbs()`, with named trait dimensions. Correlations are returned
  for the whole matrix rather than the upper triangle.

* `animal_prior_draws()`, `plot_prior_posterior()` and `plot_forest()` no
  longer take `unstructured_E`, which they accepted but never used. The
  residual correlation structure cannot change these draws: heritability uses
  only the diagonal of the residual covariance and the genetic correlations
  come from `G`. `?animal_prior_draws` says so.

* R (>= 4.0.0) is required. The package declared 3.5.0, which was wrong:
  `tools::R_user_dir()` and `cmdstanr` need 4.0.0.

## New features

* `ratones_morpho` and `ratones_pedigree` ship with the package: 35 cranial
  distances for 329 mice from the selection experiment of Penna et al. (2017),
  and the pedigree of the 525 individuals of the experiment. `data-raw/ratones.R`
  rebuilds both from the `diogro/ratones` package.
* `plot_G_ellipses()` gains an `X` argument. With `X`, the traits are replaced by
  their residuals after the fixed effects, so the scatter and the dashed
  reference ellipse have roughly the covariance `P = G + E` instead of including
  sex, age or group differences.
* `to_draws_list()` converts the arrays returned by `back_transform()` into one
  `posterior::draws_matrix` per quantity, for the `posterior` and `bayesplot`
  tools.
* `makeA()` builds the additive relationship matrix from a pedigree with the
  tabular method, accounting for inbreeding and accepting the pedigree in any
  order. `animal_design()` uses it when it is given `pedigree`.

## Improvements

* `cmdstanr` and `rstan` are suggested rather than imported, so the package
  installs without a Stan engine. A fit uses `cmdstanr` when it is available,
  falls back to `rstan` with a warning when it is not, and stops with the
  installation command when neither is present. The fallback warning used to
  suggest a malformed `install_github()` call.
* The package no longer needs `nadiv`, which was only used to build `A` from a
  pedigree: `makeA()` replaces it, so nothing has to be installed for the
  pedigree path and it works on the same R versions as the rest of the package.
  `nadiv` stays in Suggests for one test, which checks `makeA()` against it.
* The README and `?lmm_animal` state the cost that grows with the number of
  individuals, which the trait-oriented cost notes did not mention: `A` is
  dense, so each chain eigendecomposes it before sampling. The package targets
  many traits on hundreds to a few thousand individuals, and cohorts of tens of
  thousands are pointed to engines that use the sparse inverse of `A`, such as
  ASReml.
* The `?size_shape` documentation explains what the reparametrization does to
  the LKJ prior, and that a larger `lkj_prior` is often reasonable for shape
  axes. `back_transform()` and `lmm_animal_gibbs()` are documented as sharing
  the same draw layout.
* The README is reordered around the main path: quick start, usage, univariate
  and multivariate models, output, then how to build `A`, then the details.

## Bug fixes

* The `control` list now reaches both engines. With `cmdstanr`, only
  `adapt_delta` and `max_treedepth` were forwarded and every other entry was
  silently dropped; all documented sampler controls are now passed on, with
  `step_size` translated to `rstan`'s `stepsize`. Names that neither engine
  supports are an error instead of being ignored.
* The mean of the response is now restored for any design whose columns can
  represent a constant, not only when there is a literal intercept column. A
  formula such as `y ~ 0 + factor(g)` used to return coefficients on the
  centred response, missing the overall mean; they are now the group means.
  A design that cannot represent a constant (`y ~ 0 + x`) is rejected with a
  message, instead of silently returning coefficients that no longer match the
  data.
* A response trait that does not vary is now rejected with a clear message
  before the model is fitted. The univariate and multivariate Stan models
  standardize the responses by their standard deviation, so a constant trait
  used to fail inside Stan, after compiling and starting the chains.
* `animal_prior_draws()` validates `lkj_prior`, `n_traits` and `n`. A shape of 0
  returned correlations of exactly +/-1 (and a negative or non-finite shape
  returned `NaN`s with a warning), while a fractional `n_traits` or `n` was
  silently truncated and `0` failed with an unrelated error inside `matrix()`.
  `lmm_animal()` and `lm_phen_cov()` reject an invalid `lkj_prior` too, instead
  of losing every chain inside Stan.
* A count above `.Machine$integer.max` is rejected rather than coerced to `NA`
  with a warning, which used to reach the callers as an invalid matrix size.
* `plot_forest()` checks `prob` and `plot_G_ellipses()` checks `prob`, `level`
  and `n_show`, so a level of 2 or a negative number of draws is a message
  naming the argument instead of a failure inside the interval or ellipse
  calculation. `n_show` must also be a whole number: it used to be truncated,
  so `n_show = 1.5` silently drew one ellipse.
* The size and shape example in the README built its response with a fixed
  number of shape axes, so `back_transform()` failed whenever the number of
  traits was not four.
* `back_transform()` called `stats::cov2cor()` through an unqualified name,
  which `R CMD check` reported.

## Tests and documentation

* New tests cover the bundled datasets, the `X` option of `plot_G_ellipses()`,
  `to_draws_list()`, the sampler `control` pass-through of both engines, and the
  choice between `cmdstanr` and `rstan`.
* The Stan-specific branches of `animal_parameter_draws()`,
  `plot_prior_posterior()`, `plot_forest()` and `plot_G_ellipses()` are exercised
  on a real `CmdStanMCMC` fit, not only on synthetic draws, and the cache
  directory is tested through `STANIMAL_CACHE_DIR` as well as the option.
* The code blocks of the README are now evaluated as a test, so documentation
  cannot drift from the API unnoticed (skipped outside a checkout and on CRAN).
