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
