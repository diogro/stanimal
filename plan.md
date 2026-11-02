# stanimal roadmap — 0.7.0

Working plan for the next functionality release. Current state is 0.6.0: the
modelling core (animal and phenotypic models, three `A` constructions, size/shape
rotation, factor selection tooling) is in place and `R CMD check` is clean.

The theme for 0.7.0 is **datasets that are not already tidy** and **outputs a
quantitative geneticist expects but currently cannot get**: incomplete trait
records, individual breeding values, and predictive model comparison.

Two phases follow that:

- **Phase 2 — evolutionary analysis of `G`.** Everything the package is *for*:
  evolvability and constraint, where the variation sits, and reconstruction of
  selection. It takes a covariance the model already estimated and never looks at
  genotypes, because the datasets these methods exist for — comparative samples of
  skulls, of any between-species sample — essentially never come with markers.
- **Phase 3 — genomics: the multivariate effect of loci.** A second way to arrive
  at `G`, and one that also gives per-locus effect vectors, for the occasions when
  genotypes do exist. This is
  [Melo, Marroig & Wolf (2019)](https://doi.org/10.1093/jhered/esz011) and its code
  in `../mouseGrowthQTLs/`.

The same selection-reconstruction question then has to be asked of both routes, and
that comparison is the research project attached to Phase 2.

## Phase 1 — fit better data and report more

### Priority order

| # | feature | effort | value |
|---|---|---|---|
| 1 | Missing trait values | medium | very high |
| 2 | Breeding values (BLUPs) | medium | very high |
| 3 | LOO / WAIC / posterior predictive | low | high |
| 4 | Univariate `log_lik` | trivial | medium |

Rationale for the order: 1 changes what problems the package can be pointed at,
2 fills the most surprising hole in its output, and 3–4 are cheap and low risk.
3 depends on 4 (LOO for univariate fits needs its `log_lik`), so do 4 first if 3
comes early.

### 1. Missing trait values

**Today:** `Y` must be complete. `.check_numeric_matrix(Y)` rejects `NA`, and
`animal_design(na.action = stats::na.omit)` drops incomplete rows. Dropping rows
is the wrong answer when traits are correlated — the animal model can borrow
strength across traits.

**Approach.** Group rows by which traits are observed. A complete pattern keeps
the current fast path untouched. An incomplete pattern evaluates its likelihood on
the observed subvector, using the corresponding principal submatrix: with
`P = block-cholesky(P)` in transformed data, each pattern gets one Cholesky of its
own size, and the whitening is a solve against that factor.

**Changes**
- `utils.R`: allow `NA` in `Y` for the model entry points (`allow_na = TRUE` is
  already a hook there), with a clear error when a row has no observed trait or
  when an all-missing trait appears.
- `animalModelG`, `animalModelGFactor`, `phenotypicModel`, `phenotypicModelFactor`:
  add a pattern index and the per-pattern submatrix Cholesky to `transformed data`;
  the multivariate rotation path applies to complete patterns only.
- `animal_design()`: stop dropping rows by default; return the missingness mask
  alongside `Y`.
- Downstream: `pc1_variance_share()`, `size_shape()`, `scree_eigenvalues()` and
  `plot_G_ellipses()` either handle `NA` or error clearly; `plot_G_ellipses()`
  only uses `Y` for data ellipses, so it can skip incomplete rows.

**Validation**
- Simulate with a known `G`/`E`, delete 20% of entries at random, and check that
  recovery is better than the row-deletion fit on the same data, and close to the
  complete-data fit.
- Test: one trait entirely missing in a few rows; one row with a single trait;
  patterns that differ only by which trait is absent.

**Risks.** The number of distinct patterns can explode (2^K in principle) —
cap it and merge rare patterns or error, rather than silently building hundreds of
blocks.

### 2. Breeding values (BLUPs)

**Today:** `u` is marginalized out for speed and never reported. Users get `h2`
but no individual estimates, which is conspicuous for a breeding audience.

**Step 1 (0.7.0): post-processing, empirical Bayes.** A new exported function
solves the mixed-model equations using the posterior means of `G`, `E` and `beta`
and the relationship matrix `A` the fit was given:

- rotate into the eigenbasis of `A` (already computed inside the sampler, so this
  is cheap in R), solve per-trait systems, return BLUPs on the data scale;
- optional approximate standard errors from the posterior draws of `G`/`E`;
- pure R, no Stan changes, so it ships without touching the model files.

**Step 2 (later): proper posterior.** Per-draw conditional means of `u` in
`generated quantities`, which is `O(N K^3)` per draw and stores `N * K` numbers
per draw. Only worth it behind a flag (`return_breeding_values = FALSE` by
default), and it needs a memory note in the docs.

**Validation.** BLUPs from a well-specified simulation should correlate strongly
with the simulated `u`, and the empirical-Bayes solution should be close to the
posterior draw version on a small case.

### 3. LOO / WAIC / posterior predictive

`log_lik` already exists for every multivariate and phenotypic fit, and
`examples/structure_comparison.R` uses only MSE against a simulated truth.

- `loo()` method on the fit, wrapping the `loo` package and pointing at the
  `log_lik` draws.
- One posterior-predictive helper: simulate replicate responses from the
  posterior and return a data frame or a `bayesplot`-style check.
- Documentation must say what LOO means here: individual-level predictive
  adequacy, so it compares model structures, not fixed effects.

Add `loo` to `Suggests`.

### 4. Univariate `log_lik`

`animalModelUni.stan` and `animalModelUniH2.stan` are the only models without a
pointwise `log_lik`, so LOO is impossible for the fits people run most. Trivial
generated-quantities addition; check that the `h2` parametrization reports it
through the same density.

## Phase 2 — evolutionary analysis of `G`

Everything here takes a fitted covariance and never looks at genotypes. The model
supplies `G` from whatever information exists — a pedigree, a marker-based
relationship matrix, or just the phenotypes — and these functions work on that
matrix. That is what makes them usable on comparative data, which is where the
questions actually get asked.

The whole phase rests on one primitive: projecting a direction onto the axes of
`G`. Selection gradients, response vectors, evolvability, constraint and the
alignment of locus effects are all that projection, viewed differently, so they
should share one implementation.

| # | feature | effort | value |
|---|---|---|---|
| 5 | Evolvability and constraint summaries | low | high |
| 6 | `g_directions()` — where variation and selection sit | low | high |
| 7 | `reconstruct_selection()` — `beta` from a response | medium | very high |

### 5. Evolvability and constraint summaries

Port evolqg's metrics to operate across posterior draws: evolvability,
respondability, conditional evolvability, autonomy, integration, and the
mean-standardized versions. Input is a fit or a draws array, output is a tidy data
frame with one row per draw, so the posterior is summarised the same way as the
other quantities.

Cheap (pure R), and it makes stanimal the natural Bayesian backend for evolqg
users. Decide whether to depend on evolqg or reimplement the closed forms to keep
the dependency list flat — reimplementing is a few dozen lines and avoids a
weighty dependency.

The novelty over evolqg is the input: evolqg takes a point estimate of `G`, and
here every metric becomes a distribution that carries the uncertainty of the fit.

### 6. `g_directions()` — where variation and selection sit

Eigen-decompose `G` and report, per axis: the eigenvalue and its share of the
total, the projection of a supplied direction onto it (observed divergence, a
selection gradient, the vector of means, or `dz`), and the share of that direction
carried by the leading axes. Plus subspace angles between the fitted `G` and any
reference matrix, and between a direction and the constrained subspace.

This is the "size is a line of least resistance" analysis that this literature is
built on, written once and available for any fit. It is also the machinery the
scalar summaries in item 5 are built from, and the machinery that
`effect_vectors()` in Phase 3 reuses for locus effect directions.

### 7. `reconstruct_selection()` — `beta` from a response

For each posterior draw of `G`, solve `beta = G^-1 dz`, so the answer is a
distribution rather than a point: mean, credible interval, the probability that a
component has a given sign, and the cosine between the reconstructed gradient and
a reference direction (a known `beta` in a simulation, or a hypothesised direction
in an empirical analysis). `response_to_selection(fit, beta)` is the forward map,
`dz = G beta`, which is what makes the two checkable against each other.

**Numerical policy.** Inverting `G` per draw is fragile when a draw is close to
singular, and the factor models actively encourage near-singular draws. Decide and
document: a regularized inverse (eigenvalue floor, reported in the output) or
dropping the failing draws and saying how many. Never fail silently.

### The research question: which estimator reconstructs best

Reference: Marroig G, Melo DAR, Garcia G (2012). *Modularity, noise, and natural
selection*. Evolution 66(5):1506–1524.
[doi:10.1111/j.1558-5646.2011.01555.x](https://doi.org/10.1111/j.1558-5646.2011.01555.x)
(preprint: [arXiv:1112.1391](https://arxiv.org/abs/1112.1391)).

From their abstract: selection reconstruction inverts `G`, and inversion amplifies
precisely the directions that modularity and sampling or measurement noise make
least reliable; noise-controlled shrunken estimates vastly improve reconstruction
in their simulations and in a New World monkey skull analysis.

The 2019 locus-effect model changes the *estimator*, not the question, so the
study asks whether the conclusion survives the change of estimator. The hypothesis
worth testing: **the shrinkage that makes `G` estimable can make `G^-1`
untrustworthy.** Horseshoe shrinkage drives `G` towards low rank, and `G^-1` is
dominated by the smallest eigenvalues, so some of the variance that shrinkage
removes is exactly the variance the inversion divides by. Two consequences:

- there is an **optimal amount of shrinkage for reconstruction that is not the
  amount that best recovers `G`**, and the two diverge as noise grows;
- **modularity moves the optimum**. Modular `G` has small between-module
  eigenvalues, so reconstruction should over-attribute selection to between-module
  directions — the same directions where pleiotropic effects are weakest.

A negative result is informative: if the locus-based `G` reconstructs worse, the
recommendation is to report the posterior of `beta` under a weakly informative
prior, not to read a shrunken `G` as a point estimate of the covariance.

**Estimators to compare.** Sample `G` (the 2012 baseline); the 2012 noise-control
family (Ledoit–Wolf, Schäfer–Strimmer, reduced rank); pedigree or marker `G` with
the animal model; factor-structured `G`; and the locus-effect `G` from Phase 3,
with and without the polygenic term. Note the constraint that motivates the phase
split: **comparative data can only use the first four.** The mouse data can use all
of them, which is why it is the bridge between the two phases.

**Design.** Simulate from the 2019 model: modules with a controlled
within/between correlation ratio, `K` traits, `J` loci with a known effect matrix
`B` and hence a known `G`, a known `beta`, and a known `dz`. Cross a grid of
modularity (1–6 modules and the within/between ratio) × effect sparsity ×
effect-size distribution × sample size × measurement noise, with each estimator.

**Metrics.** Cosine between `beta_hat` and `beta`; per-component relative error
and sign accuracy; **alignment of the error with the smallest eigenvector of `G`**
(the mechanism, not just the symptom); coverage of the posterior intervals;
conditioning and effective rank of `G_hat`; and the angle between the reconstructed
direction and the between-module axis.

**Data.** The New World monkey skull database behind the 2012 analysis, for a
direct comparison with the published reconstructions and for the
no-genotypes case; and the mouse F3 population, where all estimators are available
at once. Both are the author's.

## Phase 3 — genomics: the multivariate effect of loci

Reference: Melo D, Marroig G, Wolf JB (2019). *Genomic Perspective on
Multivariate Variation, Pleiotropy, and Evolution*. Journal of Heredity
110(4):479–493. [doi:10.1093/jhered/esz011](https://doi.org/10.1093/jhered/esz011).
The analysis code is `../mouseGrowthQTLs/`.

The paper is not a genome scan. Its contribution is to estimate **a vector of
effects for each locus** — one number per trait, so a locus is a direction in trait
space rather than a p-value per trait — and then to characterise those vectors by
their magnitude and by their **alignment with the direction of evolutionary
divergence**. Those vectors then reconstruct ancestral trait values, decompose `G`
into per-locus contributions, and show whether the pleiotropic effects are modular.
Additive effects came out larger than dominance effects and better aligned with
selection, with the larger effects more aligned than the smaller ones.

So this phase has two jobs: give Phase 2 a second, locus-based route to `G` where
markers exist, and produce the effect vectors whose alignment Phase 2 already knows
how to measure.

### What the sibling repository contains

- `stan/SUR.stan` — seemingly unrelated regressions: K traits against J loci,
  with additive (`ad`) and dominance (`dm`) dosage vectors per individual,
  family means, effect matrices `K × J`, and `G`/`R` matrices.
- `stan/SUR_horseShoe.stan` — the same model with a horseshoe on the effects: a
  global scale `tau` per trait and a local scale `lambda` per slope.
- `stan/SUR_fHS*.stan`, `stan/finish_HorseShoe.stan` — a "finished" horseshoe
  with a slab threshold `c2`, i.e. the regularized horseshoe that stops the
  penalty from crushing large effects.
- `stan/SUR_HAL.stan` — horseshoe variant that also estimates a polygenic `G`
  and a residual `R`, so markers and pedigree are fitted together.
- `stan/mixedModelGmatrix.stan` — the pure polygenic model stanimal already has.
- `svd_mapping.R`, `genetic_effects.R`, `crossValidate.R` — effect decomposition,
  per-locus summaries and predictive checks; `R/splines.R` and
  `create_tutorial_data.R` — growth-curve prep and a crude QTL simulator.

### Proposed features, cheapest first

| # | feature | effort | value |
|---|---|---|---|
| 8 | `genotype_matrix()` and `simulate_qtl()` | low | high (foundation) |
| 9 | `locus_effects()` — joint multivariate effect matrix | high | very high |
| 10 | `effect_vectors()` — magnitude, alignment, decomposition | low | very high (needs 9) |
| 11 | `assoc(fit, M)` — fast conditional scan for screening | low-medium | medium |
| 12 | `g_contributions()` — per-locus contributions to `G` | low | high (needs 9) |
| 13 | Dominance: `genomic_relationship(type = "dominance")` | medium | medium |

**8. Foundation.** `genomic_relationship()` already filters monomorphic and rare
markers and imputes missing genotypes, but does it internally. Factor that into a
documented `genotype_matrix()` that returns aligned dosages with `Y`/`A` row
order. Add a public `simulate_qtl()` that generates genotypes from an experimental
cross, additive and dominance effects of known magnitude and direction, and a
polygenic component — the existing tests need one, and `create_tutorial_data.R` in
the sibling repo is a stopgap version.

**9. Joint multivariate effect matrix.** The paper's core model, and the `SUR`
family in the sibling repo: `K` traits against `J` loci, additive and dominance
codings, family means, and an effect matrix `B` of size `K × J` whose *columns are
the locus effect vectors*. Shrinkage on `B` with the regularized horseshoe: a
global scale per trait and a local scale per slope, which is the `SUR_fHS` variant
(global `tau[K]`, local `lambda`, slab threshold `c2`). Optional polygenic term so
the pedigree can be present at the same time, which is the `SUR_HAL` variant. Fit
per chromosome or window; document that workflow rather than pretending the whole
genome fits.

**10. Effect vectors.** The part that makes the model answer the paper's question.
For each locus, summarise its column of `B`: the norm (magnitude), the **cosine
with a supplied direction** (the divergence vector, a selection gradient, `g_max`,
or a PC axis), and the angle between loci. Across loci, eigendecompose the effect
matrix to get the dominant axes of genetic effect — the analysis in
`svd_mapping.R`. This is the same projection `g_directions()` implements for `G`, so
build it once and use it twice. Because these are posterior draws, alignment is a
distribution rather than a point estimate.

**11. Screening scan.** Still useful, but as a first pass rather than the product.
Do _not_ refit the animal model per locus — the sibling repo used `lme4qtl` for
interval mapping precisely because an MCMC scan is unaffordable. Whiten `y` and the
genotypes once with `V^{-1/2}` from the posterior means of `G` and `E`, then per
locus run a GLS across traits: `O(N)` per locus after the whitening. Marginal, and
conditional on the covariance, so document it as a screen that nominates loci, not
as an inference.

**12. Contributions to `G`.** Decompose the fitted `G` into per-locus terms from
the effect matrix (each locus contributes a rank-one piece scaled by its allele
frequency term), and compare the summed contributions with the `G` estimated from
the pedigree or the markers. This is how the paper connects the pleiotropic
distribution to the genetic covariance, and it gives a principled check that the
locus effects explain the covariance they should.

**13. Dominance.** A dominance covariance needs a second relationship matrix and
therefore a second random effect; the animal model currently has exactly one. The
marker models get dominance effects without that structural change, so do this
last, and only if the two-effect machinery is wanted for its own sake.

### Design questions for this phase

- Which shrinkage variant is the one to carry forward: plain horseshoe, the
  regularized `fHS` with its slab threshold, or `HAL` with the polygenic term? The
  paper reports both additive and dominance effects, so both codings are needed
  either way.
- Global shrinkage per trait (`tau[K]`, as in the sibling models) or one global
  scale per model? Per-trait is what a multi-trait effect matrix wants; document it
  as a prior choice the way `hs_global_scale` is documented for the factor models.
- How many loci per fit? The existing models run per chromosome. That decision
  drives memory and the recommended workflow, so it belongs in the docs.
- Ancestral reconstruction: is that a function, or an example worked through with
  the mouse data? It depends on the cross design (which lineages, which ancestral
  allele frequencies), so it may belong in `examples/` rather than in the API.
- The mouse growth, area and necropsy phenotypes plus the per-chromosome marker
  files would make excellent example data — they are yours, so confirm they can
  ship under GPL.

## Deliberately deferred

- **Repeated measures / permanent environment (`pe`)**: the biggest *modelling*
  gap, but the clean eigendecomposition trick that marginalizes `u` does not
  extend to three covariance components. Needs a design decision (nested Laplace,
  MME-based likelihood, or explicit `u` sampling) before any code.
- **Maternal effects, reaction norms, non-Gaussian traits**: all worthwhile, all
  large, and each needs its own prior discussion.
- **Sparse `A` for large pedigrees**: the dense eigendecomposition caps the model
  at a few thousand individuals. Fixing it is an architectural change, not a
  feature.
- **Factor structure for `E` alone**: `residual_structure` offers
  unstructured/diagonal/shared; a sparse `E` with unstructured `G` is a plausible
  fifth option but nobody has asked for it.
- **Phylogenetic comparative methods**: `G` here is within-population genetic
  variance, not the between-species covariance that PGLS needs. Tempting given the
  comparative datasets, but a different tool.

## Open questions

- Breeding values: empirical Bayes in 0.7.0 and posterior draws later, or wait and
  ship the full posterior once?
- Missing data: impute internally and report imputed values, or condition on the
  observed entries only and report nothing about the missing ones?
- `reconstruct_selection()` on a near-singular draw: regularized inverse with the
  floor reported, or drop the draw and report how many?
- Which direction do effect vectors align against by default: a supplied vector,
  `g_max` of the fitted `G`, or a `g_directions()` axis?
- Genomics: which shrinkage variant to carry forward, whether ancestral
  reconstruction is a function or a worked example, and whether the mouse data can
  be bundled.

## Definition of done (per release)

- `R CMD check` `Status: OK`, no new NOTEs.
- A test file per feature, with tolerances measured before they are set.
- `NEWS.md` entries under the release heading; `DESCRIPTION` version and date
  bumped.
- `AGENTS.md` updated for anything an agent would otherwise get wrong (the
  missing-data contract, what `loo()` compares, the scale of BLUPs, what an effect
  vector is, which estimator a reconstruction used).
- `README.md` gains a short section per feature, and the runnable-script index
  stays accurate — remember `test-readme.R` runs those blocks and skips in the
  tarball.
