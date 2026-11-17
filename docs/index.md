# stanimal

Bayesian animal models in [Stan](https://mc-stan.org/), for estimating
the **additive genetic covariance matrix (G)** of one or more traits
from a pedigree or a genomic relationship matrix.

stanimal is built for many traits measured on a moderate number of
individuals. Covariances, standard deviations, heritabilities and fixed
effects come back on the original scale of the data, the relationship
matrix enters through its eigendecomposition (computed once per fit),
and incomplete trait records are used as they are, without dropping
anyone.

## Install

``` r

install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev",
                                       getOption("repos")))
cmdstanr::install_cmdstan()

# install.packages("pak")
pak::pkg_install("diogro/stanimal/package")
```

Requires R 4.0.0 or later. Fits use CmdStan through cmdstanr by default,
with rstan as a fallback.

## Quick start

The package ships two real datasets measured on the same cranial
landmarks: `ratones_morpho` (329 mice from a selection experiment, with
their pedigree) and `opossum` (553 opossums, with a genomic relationship
matrix).

``` r

library(stanimal)
traits <- c("IS_PM", "NSL_NA", "PT_BA", "BA_OPI")
opossum[traits] <- log(opossum[traits])

design <- animal_design(cbind(IS_PM, NSL_NA, PT_BA, BA_OPI) ~ sex + age,
                        opossum)
fit <- lmm_animal(design$Y, design$X,
                  opossum_grm[rownames(design$Y), rownames(design$Y)])
fit$summary(c("h2", "corrG"))
draws_mean(fit, "G")
```

## What it fits

| Goal | Call |
|----|----|
| One trait, heritability | `lmm_animal(y, X, A)` |
| Several traits, covariance structure | `lmm_animal(Y, X, A)` |
| Many traits, a sparse low-rank `G` | `lmm_animal(Y, X, A, genetic_structure = "factor")` |
| `G` and `E` driven by the same latent factors | add `residual_structure = "shared_factor"` |
| Many correlated traits on one size gradient | `lmm_animal(Y, X, A, correlation_prior = "size")` |
| Phenotypic covariance only, no pedigree | `lm_phen_cov(Y, X)` |
| Missing trait values | leave the `NA`s in `Y`; `imputed_data(fit)` returns the completed data |
| Each individual’s breeding value | `breeding_values(fit)` |
| Compare two models | `loo::loo_compare(fit_loo(a), fit_loo(b))` |
| Check that a model reproduces the data | `posterior_predictive(fit)` |
| A quick approximate fit | `lmm_animal_laplace(Y, X, A)` |
| Relationship matrices | `makeA(pedigree)`, `genomic_relationship(M)` |

## Where to go next

- The **articles** work through complete analyses, comparing each
  posterior with a simulated truth or with the data: [one
  trait](https://diogro.github.io/stanimal/articles/univariate.md),
  [several
  traits](https://diogro.github.io/stanimal/articles/multivariate.md),
  [real mouse
  data](https://diogro.github.io/stanimal/articles/ratones.md), [a
  sparse factor
  `G`](https://diogro.github.io/stanimal/articles/factor-model.md),
  [choosing between covariance
  structures](https://diogro.github.io/stanimal/articles/structure-comparison.md),
  [incomplete
  records](https://diogro.github.io/stanimal/articles/imputation.md),
  [breeding
  values](https://diogro.github.io/stanimal/articles/breeding-values.md),
  [phenotypic
  covariance](https://diogro.github.io/stanimal/articles/phenotypic-covariance.md)
  and its [factor
  model](https://diogro.github.io/stanimal/articles/factor-phenotypic-covariance.md),
  [the Laplace
  approximation](https://diogro.github.io/stanimal/articles/laplace.md)
  and [the Gibbs
  sampler](https://diogro.github.io/stanimal/articles/gibbs.md).
- The **reference** documents every function and dataset.
- The [README](https://github.com/diogro/stanimal#readme) has the
  mathematics of each model, and
  [AGENTS.md](https://github.com/diogro/stanimal/blob/master/AGENTS.md)
  is a compact guide to the interface and to the prior choices the
  function signatures leave to you.

stanimal is under development, and the interface may still change.
