# stanimal roadmap — 0.7.0

Working plan for the next functionality release. Current state is 0.6.0 on
`master`, with 0.7.0 under way on `feature/0.7.0`: the
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

### 1. Missing trait values — DONE (0.7.0)

**Was:** `Y` had to be complete. Dropping rows is the wrong answer when traits
are correlated — the animal model can borrow strength across traits, and a row
deleted for one missing trait also loses the traits that *were* measured.

**What shipped.** Two mechanisms, because only one of them is valid everywhere.

- **Phenotypic models: per pattern.** No relatedness, so the records are
  independent and a row's missing values integrate out analytically: its density
  is the normal on the observed subvector, with the matching principal submatrix
  of `P`. Rows are grouped by which traits they observe; each pattern gets one
  Cholesky of its own submatrix and adds no parameters. A pattern with nothing
  observed has no likelihood and is an error.
- **Animal models: data augmentation.** This is where the plan was wrong. With
  relatedness a gap is informed by every relative's observed traits, so the
  observed-data marginal does not factorize over rows: restricted to a pattern it
  is `A[I,I] ⊗ G[O,O] + I ⊗ E[O,O]`, whose factor needs an eigendecomposition of
  `A[I,I]` per pattern — `O(N³)` each. The eigen-rotation does not escape this
  either, because a missing original value is a linear combination of *all*
  rotated records. So the gaps became parameters on the standardized response
  scale (`Y_mis`, reported), with a flat prior, which makes the posterior over
  `G`, `E` and `beta` the exact marginal over the observed data. The contribution
  of a gap is one outer product with a column of `Q`, not a re-rotation, so an
  evaluation stays `O(N·n_missing)`.

**Also changed**
- `animal_design()` no longer drops rows that miss only the response; `na.action`
  now applies to predictors, which is the only place it can.
- Univariate fits take missing records too, so the API does not change shape
  between one trait and several.
- `lm_phen_cov()` gained the pointwise `log_lik` it was documented as having.
- rstan reads a one-element R vector as a scalar, so the pattern index arrays are
  passed through `as.array()`.

**Validation, as run**
- Density identity, the decisive check: the pointwise terms at one draw equal the
  joint density computed independently in R from the same draws — per pattern on
  the observed subvector for the phenotypic model, and through
  `kronecker(G, A) + kronecker(E, I)` for the animal model. Agreement to 1e-6.
- 20% of entries deleted at random, `N = 200`, `K = 3`: `P` from the missing-data
  fit is 2.4% from the complete-data fit, while the row-deletion fit (which keeps
  only 106 of 200 rows) is 21.5% away.
- Both factor structures and the shared-factor model take gaps.

**Risk that materialised.** Not pattern explosion — that one is bounded by `N`
and only affects the phenotypic models. The risk is instead that a diagonal
residual with a *small* sample and many gaps is weakly identified: 60 individuals
with 20 gaps diverged on half the iterations, where 200 individuals with 60 gaps
is clean. Documented in `?lmm_animal` and `AGENTS.md` rather than fixed.

**Left over, now closed.** `plot_G_ellipses()` and `size_shape()` assumed a
complete `Y` and errored rather than dropping rows. Both take gaps now.

**Gaps downstream of the fit: `size_shape()`.** The two functions that earlier
comment called "modals" do consume `G` and `P`, but the size/shape
reparametrization consumes `Y` — and it is the one step that *cannot* work from
a gap, because size is a mean of the traits and shape is a difference between
them. So the response is completed before the transformation: a
maximum-likelihood multivariate normal estimated from the values that were
measured (EM, which is closed-form per pattern), and each gap replaced by its
conditional mean given the traits that individual has.

Three things worth remembering:

- The formulas are linear in the traits, so the conditional *covariance* of a
  row's gaps maps straight to the spread of its size and shape scores. That is
  where `scores_sd` comes from, and it is why the honest answer needed no
  simulation: an individual missing the trait that defines the smallest shape
  axis ends up with a score whose standard deviation is the size of that axis.
- The earlier idea of dropping incomplete rows instead was worse than it looked:
  it would have silently changed what "size" means between individuals, because
  a mean of three traits is not a mean of four.
- Point completion again, but with a sharper caveat than `imputed_data()`: here
  the *transformation* is also a point estimate, so a fit on the scores treats
  everything as data. `AGENTS.md` says so, because it is the kind of thing a
  later reader will otherwise assume was handled.

The diagnostic in `animal_design()` uses the same EM covariance and is treated as
advisory: a pattern it cannot complete is not a reason to fail a design.

### 2. Breeding values (BLUPs) — STEP 1 DONE (0.7.0)

**Today:** `u` is marginalized out for speed and never reported. Users get `h2`
but no individual estimates, which is conspicuous for a breeding audience.

**Step 1 (0.7.0): post-processing, empirical Bayes.** A new exported function
solves the mixed-model equations using the posterior means of `G`, `E` and `beta`
and the relationship matrix `A` the fit was given:

- rotate into the eigenbasis of `A` (already computed inside the sampler, so this
  is cheap in R), solve per-trait systems, return BLUPs on the data scale;
- optional approximate standard errors from the posterior draws of `G`/`E`;
- pure R, no Stan changes, so it ships without touching the model files.

**What shipped.** `breeding_values(fit, se = TRUE, ndraws = 100, seed = NULL)`
returns the conditional mean `(A ⊗ G) V^{-1} (y − Xβ)` with
`V = A ⊗ G + I ⊗ E`, at the posterior means of `G`, `E` and `beta`, as
individuals × traits on the response scale. With `A = Q diag(λ) Q'` the rotated
records are independent and record *i* has covariance `λ_i G + E`, so the `N × N`
system is never formed: it is one `k × k` solve per individual, and the
conditional variance `Var(u | y)` comes back through the same rotation. From that
fall out the two quantities a breeder actually reads:

- `se`, the square root of the conditional variance at the posterior means;
- `reliability` = `1 − Var(u_{nk} | y) / G_{kk}`, the share of a trait's additive
  variance the prediction accounts for, in [0, 1].

**Three decisions worth remembering**

- *Missing values cost nothing.* The gaps are filled with their posterior mean
  from `imputed_data()` and the prediction made from that. This is exact, not an
  approximation: the predictor is linear in the response, so `E[u | y_obs]` is
  the predictor applied to `E[y | y_obs]`. The alternative — a separate
  conditional for each pattern — buys nothing.
- *Univariate fits go through the same code.* `G` and `E` are rebuilt from
  `sigma_G` and `sigma_E` when the fit reports standard deviations instead of
  covariances, so there is one solver rather than two.
- *The step-1 shortcut is stated, not hidden.* The `se` and `reliability` are
  conditional on the posterior means, so the covariance uncertainty is missing
  from them. That is what makes step 2 worth doing rather than merely tidier.

**Validation, as run**
- The solve against the definition, written out with `kronecker()` products:
  values, conditional variance and reliability agree to 1e-10. This caught the
  one real trap, which is that `vec()` of an individuals × traits response has
  covariance `kron(G, A)` and not `kron(A, G)` — the wrong choice gives the right
  numbers in the wrong order and fails only on ordering.
- With `A = I` the prediction reduces to shrinking each residual by
  `G (G + E)^{-1}`, which is the case a reader can check by hand.
- On a simulation, the predictions correlate above 0.6 with the simulated `u`,
  and their RMSE against the truth is lower than the RMSE of each individual's
  own deviation from the fixed effects. The same test with six entries blanked
  out of a 120-individual response predicts the individuals whose trait was
  never measured, and gives them lower reliability than those that were.
- `examples/breeding_values.R` on the mouse data: reliability 0.15–0.52 across
  three traits, and only 4 of the top 10 by breeding value are in the top 10 by
  measurement. That gap is the point of the feature.

**Step 2 (later): proper posterior.** Per-draw conditional means of `u` in
`generated quantities`, which is `O(N K^3)` per draw and stores `N * K` numbers
per draw. Only worth it behind a flag (`return_breeding_values = FALSE` by
default), and it needs a memory note in the docs. The validation the plan asked
for — empirical Bayes against a per-draw posterior — needs that work or the Gibbs
sampler, which samples `u` explicitly but does not report it.

### 3. LOO / WAIC / posterior predictive — DONE (0.7.0)

**What shipped.** `fit_loo(x, criterion = c("loo", "waic"))` hands the `log_lik`
draws to the `loo` package with the chain structure kept, so
`loo::relative_eff()` has what it needs; `loo(fit)` and `loo_compare()` work on
the result. `posterior_predictive(fit, ndraws)` simulates the response from each
posterior draw and returns replicates that `print()` (observed means, SDs,
correlations against the replicates, with a tail probability) and `plot()`
(histograms with the observation marked) summarise.

**The documentation point the plan asked for, stated properly.** LOO here
compares model structures, and the reason is sharper than "it is not
individual-level adequacy": the terms come from the model block, one per record
of the *rotated* response, which is where the records are independent. Their sum
is the exact joint log density, so the criterion is valid, but with a
non-identity `A` a term is not one animal. The absolute ELPD also carries a
constant from the standardization, which cancels between fits to the same
response. Both facts are in `?fit_loo` and `AGENTS.md`.

**A correction to the plan's premise.** It said `log_lik` already existed for
every multivariate and phenotypic fit. It did not exist for the phenotypic
models: `phenotypicModel.stan` and `phenotypicModelFactor.stan` never computed
it, although `AGENTS.md` claimed they did. Fixed there first, per pattern.

**Mechanics worth remembering**
- `loo` is a Suggests, so `S3method(loo, CmdStanMCMC)` in NAMESPACE makes the
  package fail to load (`object 'loo' not found`). The methods are registered in
  `.onLoad` instead, with `setHook(packageEvent("loo", "onLoad"))` for the case
  where `loo` loads afterwards. `pkgload::load_all()` does not run `.onLoad`, so
  the test calls it explicitly.
- Both fitters now attach `stanimal_data` (the response as passed in, `X`, and
  `A`) to the fit. The posterior predictive needs it, and the breeding-value pass
  will too.
- `posterior_predictive()` slices `beta[, , j]` and `G[, , j]`; with one trait R
  drops those to vectors, so every slice goes through `matrix(...)`.

**Not done.** No `bayesplot` dependency: the plot method uses ggplot2, which was
already an Import. `examples/structure_comparison.R` still compares by MSE; it
could now use `loo_compare()` instead, which is a better demonstration than a
simulated-truth MSE and is the obvious next example to rewrite.

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
