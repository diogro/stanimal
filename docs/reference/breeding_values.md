# Estimated breeding values (BLUPs)

Predicts the additive genetic value of each individual from a fitted
animal model, on the response scale. The model itself integrates the
breeding values out and so never reports them; this recovers them by
post-processing the fit, which gives the exact conditional mean for
complete responses at posterior mean `G`, `E` and `beta`. With missing
responses, point predictions use posterior imputations and are an
approximation at those plug-in parameters.

## Usage

``` r
breeding_values(fit, se = TRUE, ndraws = 100, seed = NULL)

# S3 method for class 'stanimal_breeding_values'
print(x, n = 4L, ...)
```

## Arguments

- fit:

  A fit from
  [`lmm_animal`](https://diogro.github.io/stanimal/reference/lmm_animal.md).

- se:

  Logical: also report the conditional standard error of each prediction
  and its reliability. Default `TRUE`.

- ndraws:

  Number of completed datasets to average over when the response had
  missing values; ignored when it is complete.

- seed:

  Optional seed for those completions, so that a run is reproducible.
  The sampler's seed does not cover this step.

- x:

  An object from `breeding_values()`.

- n:

  Number of individuals to show. Default 4.

- ...:

  Unused; for method compatibility.

## Value

An object of class `stanimal_breeding_values`: a list with `values`
(individuals by traits, on the response scale), `se` and `reliability`
(the same shape, or `NULL` when `se = FALSE`), the `G` and `E` matrices
used, and the individual and trait labels.
[`print()`](https://rdrr.io/r/base/print.html) shows the predictions of
the first few individuals with their standard errors and reliabilities.

## Details

Write the model as \\y = X\beta + u + e\\ with \\u \sim N(0, A \otimes
G)\\ and \\e \sim N(0, I \otimes E)\\. The predicted breeding values are
the conditional means \$\$\hat u = (A \otimes G) V^{-1} (y - X\beta),
\quad V = A \otimes G + I \otimes E,\$\$ evaluated at the posterior
means of `G`, `E` and `beta`. This is therefore an empirical Bayes
estimate rather than a posterior: it does not propagate the uncertainty
in the covariance parameters, and the standard errors below are
conditional on them.

`se` is the square root of \\\mathrm{Var}(u \mid y\_\mathrm{observed})\\
at the same plug-in parameters. For complete responses it is based on
\\A \otimes G - (A \otimes G)V^{-1}(A \otimes G)\\. For incomplete
responses it accounts for the conditional uncertainty in the gaps rather
than treating the imputations as observed. `reliability` is \\1 -
\mathrm{Var}(u\_{nk} \mid y\_\mathrm{observed})/(A\_{nn}G\_{kk})\\, the
share of the individual's prior additive variance (which is
\\A\_{nn}G\_{kk}\\, larger than \\G\_{kk}\\ for an inbred individual)
that the prediction accounts for, between 0 and 1. A value near 1 means
the individual's relatives and its own records pin its breeding value
down; near 0 means the prediction is close to the population mean. These
are the usual BLUP quantities, so they are directly comparable with what
a classical animal-model program reports.

Predictions are deviations from the fixed effects, so they are on the
scale of the trait and can be read as "this individual is above or below
the mean of its group by this much".

Missing values are handled by averaging the completed datasets from
[`imputed_data`](https://diogro.github.io/stanimal/reference/imputed_data.md)
and predicting from that mean, which is not the approximation it looks
like: the predictor is linear in the response, so plugging in \\E\[y
\mid y\_\mathrm{observed}\]\\ gives \\E\[u \mid y\_\mathrm{observed}\]\\
at fixed `G`, `E` and `beta`. The imputed draws each use their own
posterior draw of those parameters, so the match with the plug-in values
is approximate and depends on the selected draws. The standard errors
and reliabilities use the observed-data conditional variance at the
plug-in parameters, not the variance across the imputed datasets and not
full posterior uncertainty. Computing them with \\m\\ gaps requires a
dense \\m\\ by \\m\\ matrix and its Cholesky factorization; use
`se = FALSE` to skip this cost when only predictions are needed.

Because `A` enters through its eigendecomposition, this is also a direct
check on a fit: with `A` equal to the identity the predictions are
ordinary ridge shrinkage of each individual's own residual, and with
relatedness they borrow from relatives.

## See also

[`lmm_animal`](https://diogro.github.io/stanimal/reference/lmm_animal.md),
[`imputed_data`](https://diogro.github.io/stanimal/reference/imputed_data.md),
[`draws_mean`](https://diogro.github.io/stanimal/reference/draws_mean.md)

## Author

Diogo Melo

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- lmm_animal(Y, X, A, iter = 1000, warmup = 500,
                  chains = 2, cores = 2, refresh = 0)
ebv <- breeding_values(fit)
ebv                                   # first individuals, with reliability
head(ebv$values)                      # individuals x traits
range(ebv$reliability)

# The predicted values are on the trait's scale, so they can be ranked
order(ebv$values[, 1], decreasing = TRUE)[1:5]
} # }
```
