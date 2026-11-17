# Examples

The worked examples of each model are vignettes: see the articles on the
[package website](https://diogro.github.io/stanimal/), or
`browseVignettes("stanimal")` after installing. Their sources are in
`package/vignettes/`.

This directory keeps the longer studies, which fit many models or install other
packages and are not suited to a vignette. Run them from the project root:

```sh
Rscript examples/package_comparison.R
```

| Script | Shows |
|---|---|
| `multivariate_power_simulation.R` | Parallel simulation study of genetic-correlation recovery across sample sizes, writing summaries to `results/` |
| `package_comparison.R` | The same bivariate animal model fitted with `stanimal`, `MCMCglmm`, `BGLR`, `sommer`, `lme4breeding` and ASReml-R, with estimates compared with the simulated truth, minimum ESS of `G` and timings |
| `trait_comparison.R` | The same comparison for a user-chosen number of traits (first argument, default 5), reporting relative MSE of `G`, `E`, `h2` and `corrG`, minimum ESS of `G` and timings |

All scripts need `stanimal` installed and CmdStan available
(`cmdstanr::install_cmdstan()`). The two comparison scripts install `MCMCglmm`,
`BGLR`, `sommer` and `lme4breeding` if they are missing, and check for a
licensed ASReml-R without trying to install it.
