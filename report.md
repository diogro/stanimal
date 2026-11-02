# stanimal package review

## Summary

The documentation split is a clear improvement. [`AGENTS.md`](./AGENTS.md)
now focuses on choosing and using models, while
[`DEVELOPMENT.md`](./DEVELOPMENT.md) collects repository layout, test/build
commands, and development pitfalls. The earlier review findings about custom
Stan priors, covariance-fit metadata, and scale wording have been addressed.
The main remaining issue is that one metadata promise in `AGENTS.md` is not
true for the `rstan` backend.

## Finding

### Fit trait-name metadata is absent with `rstan`

**Severity: Low**

`AGENTS.md` states that both fitters preserve `x_names` and `y_names`
([`AGENTS.md`](./AGENTS.md#L67)). In both `lmm_animal()` and `lm_phen_cov()`,
those attributes are assigned only in the `cmd_stan` branch. The `rstan`
branches rename fixed-effect draws but do not attach either attribute
([`lmm_animal.R`](./package/R/lmm_animal.R#L268),
[`lm_phen_cov.R`](./package/R/lm_phen_cov.R#L166)). As a result,
`fit_draws()` cannot recover the response's trait names from an `rstan` fit and
uses generic names such as `trait1`, even when the response columns were named.
The current `rstan` integration tests check the renamed beta variables but do
not check trait labels.

**Recommendation:** either attach `x_names` and `y_names` in the `rstan`
branches too, or document that those attributes are only kept by CmdStan fits.
Add an `rstan` test that checks `fit_draws(fit)$G` (or its univariate
equivalent) has the original trait names.

## Documentation notes

- The user guide is still fairly substantial and overlaps some of the README,
  but its model-prior and factor-identifiability advice is useful to agents.
  The new maintainer guide gives a clear home for implementation-specific
  details; consider linking to README sections for longer user explanations
  rather than repeating them if the guide grows further.
- `DEVELOPMENT.md` hardcodes `stanimal_0.6.0.tar.gz` in its build/check
  command. That works for the current `DESCRIPTION`, but the command will go
  stale after a version bump; consider updating it alongside the version or
  using a version-independent command.

## Validation

- Built the current package from `package/` and ran
  `R CMD check --no-manual --no-build-vignettes`: **Status: OK**, including
  tests and package examples.
- `git diff --check` passed.
- CI-related findings remain excluded as requested.
