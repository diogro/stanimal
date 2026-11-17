# Completed datasets from a fit with missing values

Draws the response with its gaps filled in, one completed dataset per
posterior draw. The observed entries are the data and are the same in
every draw; only the cells that were not measured vary, which is how the
uncertainty about them shows up.

## Usage

``` r
imputed_data(fit, ndraws = 100, seed = NULL)

# S3 method for class 'stanimal_imputation'
print(x, ...)
```

## Arguments

- fit:

  A fit from
  [`lmm_animal`](https://diogro.github.io/stanimal/reference/lmm_animal.md)
  or
  [`lm_phen_cov`](https://diogro.github.io/stanimal/reference/lm_phen_cov.md)
  whose response had missing values.

- ndraws:

  Number of completed datasets to draw, taken evenly across the kept
  draws.

- seed:

  Optional seed for the draws that are simulated rather than read from
  the sampler. It only matters for a phenotypic fit; an animal fit's
  gaps are already parameters.

- x:

  An object from `imputed_data()`.

- ...:

  Unused; for method compatibility.

## Value

An object of class `stanimal_imputation`: a list with `Y` (individuals
by traits by draws, with draws last, as in
[`fit_draws`](https://diogro.github.io/stanimal/reference/fit_draws.md)),
the `observed` response it came from, the `gaps` as a data frame of
individual and trait positions, and the trait names.
[`print()`](https://rdrr.io/r/base/print.html) reports the imputed value
of each gap with its 95% interval.

## Details

An animal fit's gaps are parameters, so each draw is already a draw from
the posterior of the missing values and is reported as it is. A
phenotypic fit has no such parameters – the missing values were
integrated out – so each gap is drawn from the normal conditional
distribution given that row's observed traits, using `P` and `beta` at
the same draw. The two give the same kind of object, which is why one
function covers both.

Use these for imputation-style follow-up analyses, or to see what the
model believes about the cells it never saw. Do not collapse them to a
single filled-in matrix and treat it as data: the spread across draws is
the point.

## See also

[`posterior_predictive`](https://diogro.github.io/stanimal/reference/posterior_predictive.md),
[`fit_loo`](https://diogro.github.io/stanimal/reference/fit_loo.md)

## Examples

``` r
if (FALSE) { # \dontrun{
Y[c(2, 9), 1] <- NA
fit <- lmm_animal(Y, X, A)
imp <- imputed_data(fit, ndraws = 200)
imp                                   # each gap with its interval
dim(imp$Y)                            # individuals x traits x draws
} # }
```
