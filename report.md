# Package review: stanimal 0.7.4

Review date: 2026-11-06. Scope: package functions, Stan models, documentation,
examples and tests; CI is excluded. Findings are ordered by likely impact, not
by the effort to fix them. This is a review, not a proof of correctness of every
statistical model or every possible input.

Update: all 16 findings have been fixed and verified. The original numbering
and descriptions are retained below.

## High severity

### 1. Resolved — animal posterior predictive replicates use the wrong relationship matrix

The relationship factor is now transposed before left multiplication. Fixed
draws test the factor orientation, and 20,000 simulated multivariate replicates
verify the full individual-by-trait covariance against `A ⊗ G + I ⊗ E`.
Focused predictive and imputation tests passed, including identity,
non-unit-diagonal, univariate and phenotypic cases.

[`posterior_predictive()`](package/R/posterior_predictive.R#L115) obtains
R's upper-triangular `chol(A)` but left-multiplies independent genetic draws by
it at line 126. R's factor `R` satisfies `t(R) %*% R = A`, whereas the
resulting individual covariance is `R %*% t(R)`. For
`A = matrix(c(1, .8, .8, 1), 2)`, that is `[[1.64, .48], [.48, .36]]`
instead of `[[1, .8], [.8, 1]]`. A simulation of 200,000 draws reproduced the
wrong matrix. Trait covariance is generated separately and is not the problem.
This distorts individual variances and relatedness in animal predictive checks,
even for fully observed responses.

**Recommendation:** use the lower factor `t(chol(A))` for left multiplication
and test the empirical individual-by-individual covariance for a non-diagonal
`A`, including non-unit diagonals. The existing tests check shape and aggregate
spread, not the relationship covariance.

### 2. Resolved — breeding-value uncertainty treats imputed responses as observed

Point predictions continue to use posterior-mean imputations. Standard errors
and reliabilities now condition on the observed entries at the plug-in parameter
values, including the conditional uncertainty in the gaps. Tests compare
against a direct observed-submatrix Gaussian calculation, including a real
fit, correlated traits and relatives, and a completely unobserved unrelated
individual. Focused tests and the out-of-tree package check passed.

[`breeding_values()`](package/R/breeding_values.R#L214) fills gaps with their
posterior mean, then computes the conditional variance as though the completed
response were fully observed ([lines 63–67](package/R/breeding_values.R#L63)).
Linearity can justify a plug-in *mean* under fixed parameters, but it does not
justify omitting the conditional variance of the missing response. With `A=I`,
`G=E=1`, and one completely unobserved individual, its true fixed-parameter
conditional SD is 1 and reliability is 0; the current calculation returns
`sqrt(0.5)` and 0.5. More generally, the reported SE is too small and
reliability too high for missing records.

**Recommendation:** condition the Gaussian on observed entries when computing
the variance, or include the missing-response conditional variance by the law
of total variance. Test against a direct observed-submatrix Gaussian calculation,
including an entirely missing unrelated individual and correlated relatives.
The existing missing-data tests check ranges and comparisons, not this oracle.

### 3. Resolved — LOO and WAIC for incomplete animal responses use completed-data densities

`fit_loo()` and the package's `loo()` fit methods now refuse incomplete animal
fits, including custom models. Fits or draw matrices without response metadata
warn that their likelihood cannot be verified. Complete animal fits and
incomplete phenotypic fits remain supported. Tests verify the guard, both
backend dispatch paths and the analytical observed-data marginal versus the
completed-data density. Focused tests passed; the out-of-tree package check
returned `Status: OK` with tests run separately.

Animal models fill missing values with `Y_mis` parameters and calculate
`log_lik` from the completed, rotated response, for example in
[`animalModelG.stan`](package/inst/stan/animalModelG.stan#L237) and
[`log_lik[n]`](package/inst/stan/animalModelG.stan#L271). Their sum is the
joint density of the **observed and imputed** values, not the observed-data
marginal density at that posterior draw. [`fit_loo()`](package/R/loo.R#L93)
nevertheless passes these terms to LOO/WAIC and the documentation claims the
sum is the exact observed joint likelihood. With an entirely unobserved,
unrelated row, a latent normal density is included although the observed-data
likelihood contribution is 1; dropping that term leaves its flat-prior
imputation parameter unconstrained. Rotated terms also mix missing values across
records, so they cannot simply be treated as observed-data pointwise terms.
The phenotypic models instead marginalize missing traits by pattern and are not
subject to this specific problem.

**Recommendation:** refuse or explicitly mark LOO/WAIC for incomplete animal
fits until an observed-data predictive criterion and its factorization are
implemented. Add a test against an analytically marginalized Gaussian likelihood;
the existing incomplete-data test only verifies that a result is returned.

### 4. Resolved — four-or-more-trait ellipse plots label the wrong genetic correlations

Correlation draws are now indexed with the same trait-pair table as the plot
labels. Regression tests check the rendered panel labels and HPD intervals for
every pair with two through five traits, including single-draw inputs. Focused
plot tests passed, and the out-of-tree package check returned `Status: OK`
with tests run separately.

[`plot_G_ellipses()`](package/R/plots.R#L486) collects the upper triangle in
column-major order, but looks it up using `.trait_pairs()` in row-major order.
For four traits, entries `(1,4)` and `(2,3)` are exchanged. With known
correlations `cor(A,D)=0.1` and `cor(B,C)=0.7`, the actual plot labels A–D as
0.70 and B–C as 0.10, including their intervals.

**Recommendation:** extract each correlation using the same pair-index table
used for labels. Add a four-trait test with distinct pairwise values and verify
the plotted labels, not just the ggplot class.

## Medium severity

Findings 5–11 are resolved; the descriptions below record the original review.
The fixes include matched predictive-summary masks, explicit EM failure and
corrected prior and model-comparison guidance. The incomplete-animal LOO issue
in finding 3 is now guarded against.
The source-tree test suite and focused regression tests passed, including
factor-prior simulations checked against both Stan transformations. The final
out-of-tree package check returned `Status: OK`. The README tests warned about
an unavailable PC diagnostic for fully missing rows and about ESS capping.

### 5. Resolved — LKJ warning and documentation contradict the fitted genetic prior

The unused-argument checks in [`lmm_animal()`](package/R/lmm_animal.R#L264)
and [`lmm_animal_laplace()`](package/R/laplace.R#L129) claim `lkj_prior` is
ignored whenever `E` is diagonal. However, the unstructured genetic model
always uses it for `L_Omega_G`, including with diagonal E
([variance model](package/inst/stan/animalModelG.stan#L214),
[h2 model](package/inst/stan/animalModelGH2.stan#L237)). Consequently the
warning is false and users may believe a consequential genetic prior is not
applied. The same incorrect rule appears in
[`?lmm_animal`](package/R/lmm_animal.R#L35) and [AGENTS.md](AGENTS.md#L321).
The Laplace test explicitly expects this incorrect warning.

**Recommendation:** determine LKJ use from both genetic and residual
structures; correct the help and guide, and test the full configuration table.
The ordinary factor-G/diagonal-E and shared-factor cases really do omit LKJ.

### 6. Resolved — factor prior guidance gives an incorrect expected communality

[`?lmm_animal`](package/R/lmm_animal.R#L48),
[`?lm_phen_cov`](package/R/lm_phen_cov.R#L12),
[AGENTS.md](AGENTS.md#L309), and [NEWS.md](package/NEWS.md#L58) state that
the expected prior communality is `s^2/(1+s^2)` for slab scale `s`.
The actual row-normalized quantity depends on the random normal loadings,
local/global horseshoe scales, factor count and column decay
([Stan transformation](package/inst/stan/animalModelGFactor.stan#L155)).
The formula is not its expectation. Using the implemented prior, 15,000
simulated rows with slab and global scales 1 and decay 0.8 gave mean
communalities approximately 0.20, 0.33 and 0.47 for 1, 2 and 4 factors,
respectively, not 0.5. Shared-factor loadings are not row-normalized at all.
This guidance can materially mislead prior choice.

**Recommendation:** remove the claimed formula and describe the actual
distribution and parameter dependence; demonstrate it with reproducible prior
simulation for ordinary factor G and factor P. Explain the shared-factor
representation separately.

### 7. Resolved — missingness is not applied consistently in predictive summary comparisons

[`yrep` summaries](package/R/posterior_predictive.R#L177) compute observed
means and SDs on available data but compute replicate statistics over *all*
individuals; correlations likewise use complete observed pairs but all
replicate pairs. This contradicts the instruction to compare observed entries
and can make missingness or predictor composition look like model failure.
For observed trait values `1, 2, NA, NA` and replicate values `1, 2, 20, 20`,
the current replicate mean is 10.75 rather than the matched-mask value 1.5.

**Recommendation:** retain full `yrep` for inspection, but apply each trait's
observed mask and each pair's joint observed mask before computing statistics.
Test a non-random missingness pattern with fixed, distinguishable replicate
values.

### 8. Resolved — the size/shape EM can silently return unconverged estimates

[`.normal_em()`](package/R/size_shape.R#L67) initializes from a pairwise
covariance, which need not be positive definite, and returns after 500
iterations even when the convergence condition never held
([line 100](package/R/size_shape.R#L100)). Its callers do not inspect the
iteration count. A fixed-seed 10-by-3 incomplete dataset passed the pattern
check yet reached 500 iterations without convergence; extending the same run
to 5,000 still did not converge and changed a covariance entry by up to 3.75.
`size_shape()` nevertheless returned scores without warning; the same helper
underlies `pc1_variance_share()`. Subsequent conditional variances can be
misleading, including when negative values are clipped to zero.

**Recommendation:** use a valid positive-definite initialization, detect
non-convergence and invalid covariance estimates, and explicitly warn or fail.
Test sparse and nearly singular missingness patterns, not just well-behaved
simulations. A nonconvergent case should not be presented as an ML estimate.

### 9. Resolved — univariate breeding values use squared mean SD, not mean variance

For univariate fits, [`breeding_values()`](package/R/breeding_values.R#L192)
computes `mean(sigma_G)^2` and `mean(sigma_E)^2`; its documented plug-in is
the posterior means of `G` and `E`, which are `mean(sigma_G^2)` and
`mean(sigma_E^2)`. These differ whenever the SD posterior has nonzero variance,
and can produce different predictions, SEs and reliabilities.

**Recommendation:** square each SD draw before averaging and test with
hand-constructed draws whose mean squared SD differs visibly from squared mean.

### 10. Resolved — model-comparison documentation recommends different trait sets

[`fit_loo()` help](package/R/loo.R#L74) and
[README.md](README.md#L229) list different trait sets among suitable model
comparisons, while also correctly noting the standardization constant cancels
for the *same* response. Different response dimensions or measured values
have different likelihood targets and Jacobians; their ELPDs are not directly
comparable as written. For animal models, the rotated predictive units also
depend on the relationship matrix and eigenbasis.

**Recommendation:** restrict comparisons to the same observed response,
missingness pattern and predictive units, and for rotated animal terms the same
relationship/eigenbasis. Different fixed effects on the same response can be
compared subject to the model's other LOO assumptions.

### 11. Resolved — LKJ and posterior predictive probability explanations are misleading

[`?lmm_animal`](package/R/lmm_animal.R#L35) says LKJ(1), uniform over
correlation matrices, concentrates individual correlations near +/-1 at high
dimension. Its pairwise marginals instead concentrate more strongly around
zero as dimension increases: simulations from the package's LKJ generator
gave SDs about 0.58 at K=2 and 0.30 at K=10, and `P(|r| > .8)` about 0.207
and 0.0013. Uniformity over matrices is not uniformity of each entry.

[`posterior_predictive()` help and print output](package/R/posterior_predictive.R#L47)
also claim posterior predictive tail probabilities are roughly uniform under a
correct model. Reusing the data to fit and check a model generally makes these
probabilities conservative rather than uniformly calibrated. Extreme values
can be useful diagnostics, but the stated calibration is not generally valid.

**Recommendation:** correct both claims in help, README and printed guidance;
distinguish prior over matrices from marginal pairwise priors and posterior
predictive tail areas from uniform frequentist p-values.

## Low severity

Findings 12–16 were resolved: shape axes are checked against numerical rank,
factor loadings are relabeled, the draws-conversion example is corrected,
fully missing markers are excluded from the GRM, and the opossum and example
documentation has been corrected.

The targeted regression and README tests passed, with the existing README
ESS-capping warning. An out-of-tree build and `R CMD check --no-manual
--no-build-vignettes` returned `Status: OK`.

## Validation and limitations

The following describes validation at the time of the original review.

* The package was built out of tree and `R CMD check --no-manual
  --no-build-vignettes` returned `Status: OK`, including its test phase.
  Backend tests use `skip_on_cran()` and the README is outside the built
  tarball, so this check alone does not establish that all backend tests or
  README examples ran.
* The full source-tree `testthat::test_local()` suite passed with CmdStan and
  rstan tests and the README examples. It reported one warning: the README
  example capped ESS to avoid unstable estimates. Stan also printed
  informational messages for rejected proposals.
* The covariance, ellipse labels, missingness-mask summaries, prior-distribution
  claims, missing-record breeding-value variance, EM behavior, shape rank
  failure, relabeling and GRM edge case were checked with small independent
  calculations or executable probes.
* Standalone example scripts were reviewed but not all executed as part of
  this review. A passing package check is not evidence that the identified
  statistical invariants hold.
