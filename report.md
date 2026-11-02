# stanimal package review

## Summary

The package is broadly well structured and has substantial coverage: the built
source tarball passed `R CMD check --no-manual --no-build-vignettes`, including
the testthat suite and package examples, on this machine with CmdStan and both
Stan backends installed. The findings below are limited to issues that can
mislead users or leave documentation and regression coverage out of sync.

## Findings

### 1. Prior-versus-posterior plots can compare against the wrong prior

**Severity: Medium**

[`animal_prior_draws()`](./package/R/plots.R#L26) always generates the
unstructured-model prior: independent half-normal genetic and residual scales
and an LKJ prior for genetic correlations. When no explicit prior data frame is
provided, [`plot_prior_posterior()`](./package/R/plots.R#L89) and
`plot_forest()` reach [`.prior_posterior()`](./package/R/plots.R#L135), which
generates that prior using only the number of traits and `lkj_prior`; it does
not determine the model configuration or priors used by the fit.

This is not the prior used by all accepted fits. The univariate `"h2"`
parametrization assigns a beta prior directly to h2
([`animalModelUniH2.stan`](./package/inst/stan/animalModelUniH2.stan#L88)),
and the shared-factor model has beta priors on factor heritabilities plus
factor-loading priors
([`animalModelGSharedFactor.stan`](./package/inst/stan/animalModelGSharedFactor.stan#L175)).
Factor models also do not give `G` the unstructured LKJ prior. Even for a
standard unstructured fit, a non-default `lkj_prior` must be supplied again to
the plotting helper. The resulting plot can therefore make the data appear
more or less informative than they are relative to the actual model prior.

**Recommendation:** record the relevant fit configuration and generate a
matching prior, or reject/require explicit `prior` input for fit
configurations the helper cannot represent. Add tests for `"h2"` and factor
fits so this cannot silently regress.

### 2. The README's runnable-script inventory is incomplete

**Severity: Low**

The README calls its list the “Complete runnable scripts” inventory but omits
`package_comparison.R` and `trait_comparison.R`
([`README.md`](./README.md#L376)). Both are documented in the examples index
([`examples/README.md`](./examples/README.md#L20)). Readers who use the main
README as the package index will not discover those scripts.

**Recommendation:** add both scripts to the main README list, or change the
wording so the list is explicitly illustrative and direct readers to the full
examples index.

### 3. README code examples are skipped by the built-package check

**Severity: Low**

There is a test intended to execute README R blocks, but it explicitly skips
when the README is unavailable
([`test-readme.R`](./package/tests/testthat/test-readme.R#L88)). The README is
outside the package source directory, so it is absent from the built tarball;
consequently the successful `R CMD check` run does not execute those snippets.
The test can run from a repository checkout, but this review did not establish
that a CI job runs it there.

**Recommendation:** add a repository-checkout test step to CI, or move the
canonical executable examples into a location included in the package build
and link the README to those examples.

### 4. No in-repository CI workflow is checked in

**Severity: Low**

No `.github/workflows` files were present in the tracked repository. The test
suite is substantial and passed locally, but this checkout contains no
automated workflow demonstrating that package checks, Stan integration tests,
or the repository-only README test run on future changes.

**Recommendation:** add a CI workflow that builds and checks the package, and
run the README integration test from a repository checkout with CmdStan
available (or report its intentional skip explicitly).

## Validation notes

- Built `stanimal_0.6.0.tar.gz` from `package/` and ran
  `R CMD check --no-manual --no-build-vignettes` on it: **Status: OK**.
- A direct `R CMD check package/` attempt stopped before checks because the
  source directory does not yet contain the generated `Author` field. The
  normal build-then-check path succeeded; this is a command-path caveat, not a
  failure of the built release package.
- The repository was not modified by the review aside from this report.
