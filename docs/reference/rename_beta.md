# Posterior draws with the fixed effects named

Extracts draws from a `lmm_animal` fit made with `cmd_stan = TRUE`,
renaming the `beta` variables after the columns of `X`: named
`beta[name]` for one trait and `beta[trait,name]` for several, instead
of `beta[j]` and `beta[k,j]`. Every other variable is returned
unchanged, and if `X` had no column names the draws are returned as
`fit$draws()` gives them.

## Usage

``` r
rename_beta(fit, variables = NULL, ...)
```

## Arguments

- fit:

  A `CmdStanMCMC` object returned by `lmm_animal`.

- variables:

  Optional character vector of variables to extract, using Stan's names
  (e.g. `"beta"`, `"h2"`). Default is all variables.

- ...:

  Further arguments passed to `fit$draws()`, e.g. `format`.

## Value

A `posterior` draws object, usable with bayesplot.
