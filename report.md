# stanimal package review

## Summary

The package remains well structured, with extensive tests. The new
[`AGENTS.md`](./AGENTS.md) is useful: its model-selection table, factor
identification guidance, and repository-specific test/build notes give coding
agents practical context that is difficult to infer from function signatures.
The prior-comparison behavior flagged in the previous review is now addressed
for supported built-in models: the fit records its settings, univariate
`"h2"` priors are reproduced, and factor fits require explicit prior draws.
The README's example inventory has also been completed.

## Findings

### 1. Prior comparisons can still be wrong for custom Stan models

**Severity: Medium**

`lmm_animal()` accepts `model_file` as a replacement Stan program
([`lmm_animal.R`](./package/R/lmm_animal.R#L41)), and the loader uses that
source in place of the built-in model
([`lmm_animal.R`](./package/R/lmm_animal.R#L531)). The fit nevertheless records
the ordinary model settings, and `.fit_prior_config()` uses those settings to
generate a built-in prior
([`lmm_animal.R`](./package/R/lmm_animal.R#L343),
[`plots.R`](./package/R/plots.R#L185)). If a custom model changes the prior,
`plot_prior_posterior(fit)` may silently compare the posterior to the wrong
prior.

**Recommendation:** record that a custom model was used and require an
explicit `prior` for such fits, or otherwise make the custom model's prior
metadata explicit. Add a regression test for this behavior.

### 2. `AGENTS.md` overstates configuration metadata available on both fitters

**Severity: Low**

The guide says both model functions attach
`genetic_structure`, `residual_structure`, `parametrization`, `lkj_prior`,
`beta_prior` and `n_factors`
([`AGENTS.md`](./AGENTS.md#L56)). That metadata is attached by
`lmm_animal()`; `lm_phen_cov()` only adds `x_names` and `y_names` in its
CmdStan path and does not attach covariance configuration
([`lmm_animal.R`](./package/R/lmm_animal.R#L339),
[`lm_phen_cov.R`](./package/R/lm_phen_cov.R#L166)). This could lead an agent to
rely on attributes that are absent, or to assume phenotypic factor fits retain
their configuration.

**Recommendation:** scope the listed attributes to `lmm_animal()`, and
document separately which metadata `lm_phen_cov()` currently preserves.

### 3. The guide's opening scale claim conflicts with its factor-model detail

**Severity: Low**

The guide opens with “Everything is returned on the original data scale”
([`AGENTS.md`](./AGENTS.md#L3)), but the ordinary factor model also reports
`s_G` on the standardized response scale, as the guide later correctly notes
([`AGENTS.md`](./AGENTS.md#L174)). This could be simplified to say that
covariances, standard deviations, heritabilities and fixed effects are returned
on the original scale, while factor-specific quantities such as `s_G` and
loadings have their documented model scales.

## Guide design note

The agent guide is valuable, especially its non-obvious factor-identification
and test/build advice. At 236 lines, however, it is long for a root-level
always-on instruction file and repeats some README material. Consider keeping
the short decision table and high-impact agent gotchas there, with links to the
README for longer explanations and examples.

## Validation

- Focused `prior-posterior` test file passed.
- Built the current package from `package/` and ran
  `R CMD check --no-manual --no-build-vignettes`: **Status: OK**, including
  tests and package examples.
- This report omits CI-related findings as requested.
