<p align="left">
  <img src="package/man/figures/logo.png" width="200" alt="stanimal logo: a lamb with a halo inside a hexagon (saint animal)">
</p>

# Stanimal

Bayesian animal models in Stan, with a focus on estimating the **additive genetic covariance matrix (G)** of one or more traits from a known relationship matrix. The intended niche is many traits measured on a moderate number of individuals: the cost that grows with the number of traits is the likelihood, and the one that grows with the number of individuals is the dense relationship matrix (see the cost notes under [Multivariate model](#multivariate-model)).

Fits use [CmdStan](https://mc-stan.org/cmdstanr/) through `cmdstanr` (default) or `rstan`.

# Install

```r
if (!requireNamespace("pak", quietly = TRUE))
    install.packages("pak")
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev",
                                       getOption("repos")))
cmdstanr::install_cmdstan()
pak::pkg_install("diogro/stanimal/package")
```

Requires R 4.0.0 or later.

# Quick start with the bundled data

The package ships the measurements and the pedigree of the mouse artificial-selection experiment of [Penna et al. (2017)](https://doi.org/10.1111/evo.13304):

- `ratones_morpho`: 329 mice with 35 cranial distances plus sex, age, generation and selection line, one row per individual;
- `ratones_pedigree`: the pedigree of the 525 individuals of the experiment, ancestors included.

The rows of `ratones_morpho` are named with the individual IDs, so a relationship matrix built from the pedigree is matched to the records automatically:

```r
library(stanimal)
traits <- c("IS_PM", "NSL_NA", "PT_BA", "BA_OPI")
ratones_morpho[traits] <- log(ratones_morpho[traits])   # morphometrics on the log scale

# SEX + AGE + LIN: the lines were selected for body weight, so their means differ
design <- animal_design(cbind(IS_PM, NSL_NA, PT_BA, BA_OPI) ~ SEX + AGE + LIN,
                        ratones_morpho, pedigree = ratones_pedigree)
fit <- lmm_animal(design$Y, design$X, design$A)
fit$summary(c("h2", "corrG"))
```

`examples/ratones.R` fits this model end to end and plots the results.

# Usage

```r
library(stanimal)

# Y: trait matrix (individuals in rows), one record per individual
# X: fixed effects model matrix (full column rank, and able to represent a
#    constant: include an intercept column). Optional: leave it out for an
#    intercept-only fit
# A: relationship matrix, same individual order as Y and X
A <- makeA(pedigree)                   # or genomic_relationship(genotypes)
fit <- lmm_animal(Y, X, A)
```

`Y` can be a vector (univariate) or a matrix (multivariate). `X` is optional: `lmm_animal(Y, A = A)` fits an intercept-only model, with the intercept named `(Intercept)` in the `beta` draws. `A` must be symmetric and positive definite; if it is only positive semidefinite (duplicated individuals, for example), fit with `blend_A = TRUE`. When `A` and the rows of `Y` are both named, `A` is put in the row order of `Y`; names that disagree are an error. Without names, the rows are assumed to be in the same order and only the dimensions are checked.

## Univariate models

Two parametrizations are available through `parametrization`:

- `"variance"` (default): estimates the genetic and residual standard deviations with half-normal priors and derives h2.
- `"h2"`: estimates a total standard deviation and h2, with a `beta(beta_prior)` prior on h2. The default `beta_prior = c(1, 1)` is flat on h2.

If the `"variance"` version has chains with low E-BFMI, `lmm_animal()` warns and suggests switching to `parametrization = "h2"`. This tends to happen when the genetic variance is weakly identified, for example with few families.

```r
fit <- lmm_animal(y, X, A, parametrization = "h2")
```

## Multivariate model

The multivariate model estimates `G` and the residual covariance `E` directly. Each is a vector of scales times an LKJ-distributed correlation matrix, with half-normal priors on the scales (of the standardized responses). There is no prior on h2; `h2 = diag(G) / diag(P)` is derived.

- `lkj_prior` sets the LKJ shape for both correlation matrices (default 2; 1 is uniform, larger shrinks toward zero correlation). A shape above 1 is deliberately recommended for many traits: at high dimensionality, the uniform prior puts most of its mass near ±1 correlations and overestimates their magnitude, which is a well-known problem when estimating genetic correlations in quantitative genetics.
- `residual_structure = "unstructured"` (default) fits a full residual covariance. `"diagonal"` assumes independent residuals, which pushes environmental covariance between traits into `G` and inflates heritability, so it is not recommended when estimating `G`.
- Breeding values are marginalized using the eigendecomposition of `A`, so the model has no per-individual random effects. This is much faster and better mixing than sampling breeding values, and the likelihood only needs the eigenvalues of `A`, which is why a singular `A` works once it is regularized with `blend_A = TRUE`.
- The per-individual likelihoods are also whitened by the residual Cholesky factor and rotated into a single eigendecomposition of $E^{-1/2}GE^{-1/2}$, so an evaluation costs $O(NK^2 + K^3)$ rather than $O(NK^3)$. The model scales to dozens of traits (e.g. 40) without per-individual Cholesky factorizations. When the eigenvalues of $E^{-1/2}GE^{-1/2}$ are nearly coincident (at initialization, or when `G` is roughly proportional to `E`), the model automatically falls back to a per-individual Cholesky evaluation, whose gradient is well behaved.
- Fitting many traits is still computationally heavy even though the likelihood is cheap: an unstructured `G` at `K = 40` has roughly 1,700 parameters, and NUTS needs many leapfrog steps, so one fit can take an hour or more. Fewer `iter`, a larger `lkj_prior`, or `residual_structure = "diagonal"` make fits faster.
- The side that does not scale is the number of individuals. `A` is dense, so the models receive all $N \times N$ of it and each chain eigendecomposes it before sampling: $O(N^3)$ time and $O(N^2)$ memory per chain, about two seconds at 1,000 individuals, half a minute at 2,000 and a minute at 4,000. That is a one-off cost, and afterwards an evaluation is only linear in $N$: with many traits, sampling dominates and the decomposition is a small fraction of the fit (at 4,000 individuals and `K = 3`, about a minute out of ten), while with few traits and short runs it can be most of it (at 4,000 individuals univariate, about 70% of 500 iterations). Either way the package targets many traits on hundreds to a few thousand individuals; for cohorts in the tens of thousands, where the dense relationship matrix is the bottleneck, an engine that works with the sparse inverse of `A` (such as [ASReml](https://asreml.locustware.com/)) is the better tool.

### Sparse factor models for G and E

When the number of traits is large, an unstructured `G` and `E` need $O(K^2)$ parameters and become expensive to estimate. Two factor-structured alternatives are available.

**Factor `G`.** `genetic_structure = "factor"` models the genetic covariance as

$$
G = \mathrm{diag}(s_G)(\Lambda\Lambda^\top + \Psi)\mathrm{diag}(s_G),
$$

with sparse loadings $\Lambda$ on the correlation scale and diagonal uniqueness $\Psi$. The loadings have regularized-horseshoe priors, and `column_decay` can favour factors in decreasing order. `n_factors` sets the number of factors; if omitted, a sublinear default is used.

```r
fit <- lmm_animal(log_traits, X, A,
                  genetic_structure = "factor", n_factors = 3)
fit$summary(c("loadings", "factor_magnitude"))
```

**Shared factors for `G` and `E`.** `residual_structure = "shared_factor"` uses the same latent factors for genetic and residual covariance:

$$
G = \Lambda\mathrm{diag}(h^2_{\text{factor}})\Lambda^\top + \mathrm{diag}(\Psi_G),
\qquad
E = \Lambda\mathrm{diag}(1-h^2_{\text{factor}})\Lambda^\top + \mathrm{diag}(\Psi_E),
$$

so each factor has a heritability. This follows the sparse factor model of [Runcie and Mukherjee (2013)](https://doi.org/10.1534/genetics.113.151217) (“Bayesian Sparse Factor Analysis of Genetic Covariance Matrices”, *Genetics*; [arXiv:1211.3706](https://arxiv.org/abs/1211.3706)). It requires `genetic_structure = "factor"`:

```r
fit <- lmm_animal(log_traits, X, A,
                  genetic_structure = "factor",
                  residual_structure = "shared_factor",
                  n_factors = 3)
fit$summary(c("factor_heritability", "G", "E"))
```

Factor fits report `loadings`, `factor_magnitude`, `factor_heritability`, `uniqueness` and `communality`. The last two are the less familiar quantities: `communality[k]` is the variance that trait `k` shares through the factors, and `uniqueness[k]` is its trait-specific remainder. In the ordinary factor model they are on the correlation scale, so `communality + uniqueness = 1`; multiplying them by `s_G[k]^2` gives the common and unique genetic variances. In shared-factor fits, `uniqueness = genetic_uniqueness + residual_uniqueness` and `factor_heritability[m]` is the genetic share of factor `m`. `factor_magnitude[m]` is the squared norm of loading column `m` across traits, not a per-trait quantity. Loadings and factor magnitudes both depend on how the factors are oriented, and nothing in the model fixes that orientation except the prior: two solutions with the same `G` can differ by a rotation, a sign flip or a permutation of the columns. So interpret `G`, `E`, the communalities and the uniquenesses, and treat a poor `Rhat` on `loadings` or `factor_magnitude` as an orientation that moved rather than as a failed fit. Lowering `column_decay` orders the columns and makes the magnitudes comparable. The factor prior reduces the number of covariance parameters, but the current Stan models still use a $K \times K$ eigendecomposition, so the factor structures help parameter count and parsimony more than the linear-algebra cost. `examples/factor_model.R` and `examples/structure_comparison.R` show both alternatives.

### Choosing the number of factors

`n_factors` is the argument to think about: too few factors underfit the covariance, and too many split one direction across several columns, which is what makes the loadings rotate. Two helpers look at the data instead of guessing:

```r
set.seed(1)
suggested <- suggest_n_factors(log_traits, X)      # parallel analysis
scree <- scree_eigenvalues(log_traits, X)          # the table behind the plot
plot(scree, reference = suggested)                 # scree plot with the cutoffs
```

`suggest_n_factors()` runs a parallel analysis (Horn 1965). It compares the eigenvalues of the residual trait correlation matrix with those of data that has no shared structure, after removing the same fixed effects, and counts the components that clear the cutoff. Because the cutoff accounts for the sample size, the number of traits and the covariates, it is more reliable than "keep the eigenvalues above one" or eyeballing the elbow.

Treat the result as an **upper bound**, not a target. On the mouse morphometric data it keeps seven dimensions, but the first two carry 28% and 11% of the variance while the seventh clears its cutoff by 0.003 — a statistic, not a factor. `attr(suggested, "share")` reports how much of the variance the retained components take, and `attr(suggested, "retained")` is the raw count, which is `0` when nothing stands above the noise at all. Read the scree plot and stop where the eigenvalues flatten; that is the number of factors you can interpret. If you only want a good `P` or `G`, a few more is harmless provided `column_decay` keeps the columns ordered.

That is what `examples/factor_phenotypic_covariance.R` works through. The scree on the 35 log traits shows two components at 28% and 11% and then a flat tail. A three-factor fit gives bimodal communalities (`Rhat` around 1.5 across four chains) — the model saying the third direction is not there. Two factors converge, but only with a decay strong enough to separate the two columns: at `column_decay = 0.8` their magnitudes were 6.0 and 6.7 with `Rhat` 1.5, while `0.5` gave 7.8 and 5.0 with `Rhat` 1.0.

`column_decay` (default `0.8`) multiplies the prior standard deviation of column `m` by `column_decay^(m - 1)`, which is what stops the factor labels from swapping between chains. With `column_decay = 1` the columns are exchangeable a priori, and `factor_magnitude` can report large `Rhat` values on a perfectly healthy fit. Falling back to the default `n_factors` warns with the number it used, because that fallback only sees the number of traits.

How strong the decay needs to be depends on how many columns have to be told apart. Few factors need *more* shrinkage than many: with two columns of comparable size, `0.8` is often too mild and their labels still swap, while a broader model wants the mild default, because `0.5^(m - 1)` shrinks the prior of the fifth column to a sixteenth of the first. So lower `column_decay` when a fit with two or three factors has `factor_magnitude` `Rhat` values well above one, and leave it at the default when there are many.

## Output

With the default `cmd_stan = TRUE`, `lmm_animal()` returns a `CmdStanMCMC` object, so the usual cmdstanr, posterior and bayesplot tools work directly:

```r
fit <- lmm_animal(Y, X, A)          # multivariate, so the fit holds G, E and h2
fit$summary(c("G", "E", "h2"))
fit$draws("G", format = "draws_matrix")

bayesplot::mcmc_trace(fit$draws("h2"))
bayesplot::mcmc_intervals(fit$draws("corrG"))
bayesplot::mcmc_nuts_energy(bayesplot::nuts_params(fit))
```

Variables are reported on the **original scale** of the data:

| Model | Variables |
|---|---|
| Multivariate | `G`, `E`, `P = G + E`, `corrG`, `h2`, `sigma_G`, `sigma_E`, `beta[trait, j]` |
| Univariate | `h2`, `sigma_G`, `sigma_E` (standard deviations), `beta[j]` |

Fixed effects are `beta[j]` (univariate) or `beta[trait, j]` for column `j` of `X`. To get draws named after the columns of `X` (for example `beta[2,age]`), use:

```r
rename_beta(fit, "beta")
```

With `cmd_stan = FALSE` the fit is a `stanfit` object from `rstan`, with `beta` already named after the columns of `X`.

### Covariance matrices as matrices

`fit$draws()` reports `G` one element per variable (`G[1,1]`, `G[2,1]`, ...), so getting a matrix back means matching and reshaping names by hand. `fit_draws()` does that once and returns arrays indexed by trait, with draws last and every dimension named:

```r
d <- fit_draws(fit)
d$G                                # traits x traits x draws
d$E
d$corrE                            # derived from the E draws
d$h2                               # traits x draws
d$beta                             # traits x terms x draws

round(draws_mean(fit, "G"), 2)     # posterior mean G, as a matrix
round(draws_mean(fit, "corrE"), 2)
draws_quantile(fit, "corrG", probs = c(0.05, 0.95))   # 90% interval per element
```

`fit_draws()` also accepts `lmm_animal_gibbs()` fits (returning its `R` as `E`), `lm_phen_cov()` fits, `posterior` draws objects and the output of `back_transform()`. Quantities the fit does not report, such as `corrE` or `P`, are derived from the ones it does.

## Formula interface

`animal_design()` builds `Y` and `X` from a formula and a data frame, for when you would rather write the model the usual way than assemble the matrices by hand:

```r
design <- animal_design(cbind(t1, t2) ~ age + sex, data = dat)
fit <- lmm_animal(design$Y, design$X, A)
```

- the left-hand side is a single response or a `cbind()` of several traits, and the columns of `Y` are named after them;
- `X` is the fixed-effects design matrix, so its columns match the `beta` draws reported by the fit; factors are dummy-coded and the intercept is kept unless the formula removes it, which is only allowed when the remaining columns still sum to a constant (the 0/1 indicators of a factor, for example, but not `y ~ 0 + x`), because the models centres the responses and report their mean back through the design;
- random-effect terms such as `(1 | animal)` are ignored, so `lme4`-style formulas can be reused unchanged — relatedness is supplied through `A`;
- rows with missing values are dropped (`na.action = stats::na.omit` by default), and the individuals kept, in row order, are returned in `attr(design, "individuals")`. Because `Y` keeps them as row names, a named `A` is reordered to match automatically:

```r
rownames(dat) <- dat$id
design <- animal_design(cbind(t1, t2) ~ age + sex, data = dat)
fit <- lmm_animal(design$Y, design$X, A)   # A is reordered to the rows of Y
```

With many traits, the `cbind()` response is easier to assemble from the column names:

```r
traits <- c("t1", "t2", "t3")
f <- reformulate(c("age", "sex"),
                 response = sprintf("cbind(%s)", paste(traits, collapse = ", ")))
design <- animal_design(f, data = dat)
```

The result feeds `lmm_animal()`, `lm_phen_cov()` and `lmm_animal_gibbs()` alike.

`animal_design()` can also build `A`, either from a pedigree (through `makeA()`) or from a genotype matrix (through `genomic_relationship()`), as long as `data` has the individual IDs as row names. `A` is restricted to the rows kept in `Y`, so ancestors and individuals without records are handled:

```r
rownames(dat) <- dat$id

# from a pedigree (first three columns: id, sire, dam)
design <- animal_design(cbind(t1, t2) ~ age + sex, dat, pedigree = ped)
fit <- lmm_animal(design$Y, design$X, design$A)

# or from genotypes (individuals in rows, markers in columns)
design <- animal_design(cbind(t1, t2) ~ age + sex, dat, genotypes = M, blend = 0.01)
fit <- lmm_animal(design$Y, design$X, design$A)
```

## Size and shape reparametrization (use this when one PC dominates)

**If the traits are morphometric and a single principal component explains most of their variance, fit size and shape scores instead of the raw traits.** In that situation the covariance is a long, thin ridge: size and shape are still estimable, but the raw parametrization samples poorly (divergences, low E-BFMI, many leapfrog steps, etc.) and the prior is spread over near-degenerate directions.

Check the share of variance in the first principal component first:

```r
pc1_variance_share(log_traits)   # 1 / n_traits means no dominant gradient; > 0.8 means fit size and shape
```

`size_shape()` rewrites the traits as log size plus orthogonal shape axes, and `back_transform()` maps the posterior draws back to the trait scale:

```r
ss <- size_shape(log_traits)                  # individuals in rows, traits in columns
ss_dat <- cbind(morphology, ss$scores)        # size, shape1, shape2, ...

# the response is every score, so the model traits match the ones back_transform()
# expects; build it from ss$traits instead of writing cbind(size, shape1, ...)
f <- reformulate(c("sex", "site"),
                 response = sprintf("cbind(%s)", paste(ss$traits, collapse = ", ")))
design <- animal_design(f, ss_dat)
fit <- lmm_animal(design$Y, design$X, A)

post <- back_transform(fit, ss)               # G, E, P, corrG, corrP and h2 on the log scale
draws_mean(post, "G")                         # posterior mean G
draws_mean(post, "h2")                        # posterior mean h2 per trait
```

`fit_draws()` accepts that same object if you want all of the draws in one list of arrays, and `to_draws_list()` reshapes it into one `posterior` draws matrix per quantity, for when you want the `posterior` or `bayesplot` tools on the back-transformed draws.

`animal_design()` warns automatically when the response has three or more traits and its first principal component explains 80% or more of the variance.

The LKJ prior is not invariant under this change of traits: it acts on correlations between size and shape, which are expected to be small because the shape axes are principal components, so a larger `lkj_prior` than for the same traits fitted raw is often reasonable — see `?size_shape`.

Size is `rowMeans(Y)` and the shape axes are the principal components of `Y - size`, so `Y = scores %*% t(loadings)` up to the means and covariance matrices map exactly as $B M B^\top$. Keeping fewer axes (`n_shape`) gives a lower-rank, singular covariance in the original trait space. Heritabilities and genetic correlations are recomputed after the transformation, since neither is invariant under a change of traits; `beta` stays on the model scale. `back_transform()` also accepts an `lm_phen_cov()` fit, in which case it returns `P`, `corrP` and `sigma_P`.

## Pedigree relationship matrix

`makeA()` turns a pedigree into the additive relationship matrix, using the tabular method, with parents allowed in any order and unknown parents coded `0`, `NA` or `""`:

```r
# first three columns: individual, sire, dam
A <- makeA(pedigree)
fit <- lmm_animal(Y, X, A)
```

Inbreeding is accounted for, so a pedigree with matings between relatives gives diagonal values above one. A parent that is not one of the pedigree's individuals is added as a founder, with a warning, rather than silently dropping the relationship; a pedigree where an individual is its own ancestor is an error.

### Related founders and metafounders

Founders are unrelated by default. If the base population has structure — multiple breeds or populations, or a known relationship among the founders — pass `Gamma`, the relationship matrix among founders, to `makeA()`:

```r
# a three-individual pedigree: founders 1 and 2 and their offspring 3
ped2 <- data.frame(id = c(1, 2, 3), sire = c(NA, NA, 1), dam = c(NA, NA, 2))
Gamma <- matrix(c(1, 0.25, 0.25, 1), 2,
                dimnames = list(c("1", "2"), c("1", "2")))
A <- makeA(ped2, Gamma = Gamma)
```

Founders are individuals with both parents unknown. `Gamma` can be matched by founder IDs through its dimnames, or by founder order when it has no names. This is the \(A^\Gamma\) construction for related founders and metafounders (Legarra et al. 2015): the founder block is initialized from `Gamma` instead of the identity, and the tabular recursion propagates those relationships.

For genetic groups or metafounders, assign each founder to a group and give `Gamma` one row and column per group:

```r
Gamma_groups <- matrix(c(1, 0.2, 0.2, 1), 2,
                       dimnames = list(c("A", "B"), c("A", "B")))
founders <- pedigree$id[is.na(pedigree$sire) & is.na(pedigree$dam)]
groups <- rep(c("A", "B"), length.out = length(founders))
A <- makeA(pedigree, Gamma = Gamma_groups, groups = groups)
```

`groups` needs one label per founder, either in founder order or named by founder ID.

`animal_design()` passes both arguments through when it builds `A`:

```r
design <- animal_design(cbind(t1, t2) ~ sex, dat,
                        pedigree = pedigree, Gamma = Gamma_groups,
                        groups = groups)
```

Only founders with both parents unknown receive founder relationships. An unknown parent of an individual that has one known parent is still treated as an unrelated founder; add it explicitly to the pedigree and assign a group if it belongs to the related base population. A genomic relationship matrix among founders can be used as `Gamma`, but it estimates realized rather than base-population relationships, so blending or shrinkage is advisable. See `?makeA` for details and the [Legarra et al. (2015)](https://doi.org/10.1534/genetics.115.177014) reference.

## Genomic relationship matrix

If you have marker genotypes instead of a pedigree, `genomic_relationship()` builds a VanRaden (2008) GRM to use as `A`:

```r
# M: individuals in rows, markers in columns, coded 0/1/2 (NA allowed)
A <- genomic_relationship(M, blend = 0.01)
fit <- lmm_animal(Y, X, A)
```

Markers with low minor allele frequency (`maf`), too much missing data (`max_missing`) or no variation are dropped, and missing genotypes are set to their expected value. Known base-population allele frequencies can be passed with `p`. If there are fewer markers than individuals, or duplicated genotypes, the GRM is singular; a small `blend` mixes in the identity to make it positive definite.

## Phenotypic covariance only (P)

With no pedigree or genomic data, `lm_phen_cov()` estimates just the phenotypic covariance matrix `P`. It takes the same arguments and priors as `lmm_animal()` (without `A`), and returns the same kind of object, with `P`, `corrP`, `sigma_P` and `beta`. The likelihood is a single whitened triangular solve, $O(NK^2 + K^3)$, so it also scales to dozens of traits:

```r
fit <- lm_phen_cov(Y, X, lkj_prior = 2)
fit$summary("P")
```

For many traits, the phenotypic covariance can instead use the same sparse factor structure as the animal model:

```r
fit <- lm_phen_cov(Y, X, covariance_structure = "factor", n_factors = 1)
fit$summary(c("P", "loadings", "factor_magnitude"))
```

`covariance_structure = "factor"` uses the regularized-horseshoe loadings and `column_decay`, and `n_factors` falls back to a sublinear default when it is omitted (here the two traits leave room for one factor). The factor fit reports `loadings`, `factor_magnitude`, `uniqueness` and `communality` in addition to `P`, `corrP` and `sigma_P`. `examples/factor_phenotypic_covariance.R` fits the full 35-trait mouse covariance matrix this way on the log-trait scale, and checks the fitted `P` against the OLS residual covariance.

Keep the factors on the observed trait scale. Rotating first (for instance to the size/shape scores of `size_shape()`) spreads every loading over all the traits, which removes the sparsity the prior is there to find, and turns the low-rank structure into the rotation itself: the score covariance is then almost diagonal, so the factor decomposition is left weakly identified.

## Data scaling and back-transformed coefficients

All models standardize the response and divide continuous predictors by their standard deviation (constant and 0/1 columns are not scaled) for better sampling. Coefficients are back-transformed, so `beta` is in the original units of `X` and `Y`, and the mean of the response is restored through the design: onto an intercept column, or split across the indicator columns of a factor when the intercept is dropped.

## Compilation cache

Models fit with cmdstanr are compiled once and the executable is reused in later R sessions. The cache is in `tools::R_user_dir("stanimal", "cache")` by default; change it with `options(stanimal.cache_dir = "path")` or the `STANIMAL_CACHE_DIR` environment variable. Use `stan_cache_dir()` to see the location, and run `stanimal_clear_cache()` after upgrading CmdStan.

## Other options

- `cmd_stan = FALSE` uses `rstan` instead of `cmdstanr`. Both engines are optional: without `cmdstanr` a fit falls back to `rstan` with a warning, and with neither installed it stops and says how to install one.
- `blend_A = TRUE` replaces `A` with `0.99 * A + 0.01 * I`.
- `model_file` runs a custom Stan file in place of the default model.
- Further arguments (`...`) are passed to the sampler, for example `refresh = 0`.

## Gibbs sampler

`lmm_animal_gibbs()` is a quick single-chain Gibbs sampler for simulations and sanity checks. Its priors (inverse-Wishart) differ from those in `lmm_animal()`, so the results are not directly interchangeable.

# Examples

Complete runnable scripts are in [`examples/`](examples/): `univariate.R`, `multivariate.R`, `factor_model.R` (sparse factor `G`), `structure_comparison.R` (the five covariance structures), `ratones.R` (real mouse data), `phenotypic_covariance.R`, `factor_phenotypic_covariance.R` (the full 35-trait factor P model, on the log-trait scale), `gibbs.R` and `multivariate_power_simulation.R`. They compare the posterior with a simulated truth or with the data and write figures to `examples/figures/`; see [`examples/README.md`](examples/README.md) for what each one covers.

## A complete example

```r
set.seed(1)
n_fam <- 40; fam_size <- 10; n <- n_fam * fam_size
family <- rep(seq_len(n_fam), each = fam_size)
A <- outer(family, family, "==") * 0.5   # half-sib families
diag(A) <- 1

G <- matrix(c(1.0, 0.5, 0.5, 1.5), 2)
E <- matrix(c(1.0, 0.3, 0.3, 1.0), 2)
X <- cbind(Intercept = 1, age = rnorm(n, 50, 10), sex = rbinom(n, 1, 0.5))
B <- cbind(c(10, 0.1, 1), c(-5, 0.05, -0.5))
Y <- X %*% B +
  t(chol(A)) %*% matrix(rnorm(n * 2), n) %*% chol(G) +
  matrix(rnorm(n * 2), n) %*% chol(E)
colnames(Y) <- c("t1", "t2")

fit <- lmm_animal(Y, X, A)
fit$summary(c("G", "corrG", "h2"))
rename_beta(fit, "beta")
```

# Notes

- The models assume one record per individual.
- Under development: the API may still change.
