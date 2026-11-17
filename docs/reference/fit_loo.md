# Leave-one-out cross-validation for a fit

Computes the expected log pointwise predictive density with
Pareto-smoothed importance sampling (or WAIC), leaving out one
individual at a time. Animal fits with missing trait records are
supported.

## Usage

``` r
fit_loo(
  x,
  criterion = c("loo", "waic"),
  r_eff = NULL,
  unit = c("individual", "record"),
  eigen_A = NULL,
  transform = NULL,
  ...
)

# S3 method for class 'CmdStanMCMC'
loo(x, ...)

# S3 method for class 'stanfit'
loo(x, ...)
```

## Arguments

- x:

  A fit object.

- criterion:

  `"loo"` for Pareto-smoothed importance sampling, or `"waic"` for the
  widely applicable information criterion.

- r_eff:

  Relative effective sample sizes, one per observation. Computed from
  the within-chain correlation of the pointwise likelihoods when left
  `NULL`, which is what
  [`relative_eff`](https://mc-stan.org/loo/reference/relative_eff.html)
  recommends.

- unit:

  For animal fits, what one pointwise term is. `"individual"` (default)
  leaves out one individual: the term is the density of its observed
  traits given every other individual, computed after the fit from the
  draws of `G`, `E`, `beta` and, with gaps, `Y_mis`. `"record"` uses the
  model's own `log_lik`, one term per record rotated by the eigenvectors
  of `A`, which needs a complete response. Ignored for phenotypic fits
  and draws matrices, whose terms are already one per individual or are
  taken as given.

- eigen_A:

  Optional eigendecomposition of the relationship matrix the fit used
  (`attr(fit, "stanimal_data")$A`, after any reordering and `blend_A`),
  as a list with `Q` and `lambda`, to skip the \\O(N^3)\\ decomposition
  for `unit = "individual"`. It is checked against that matrix.

- transform:

  Optional `size_shape` object, when `x` was fitted to its `scores`. The
  terms are then moved to the original traits, so the result can be
  compared with a fit to those traits (see Details).

- ...:

  Further arguments passed to
  [`loo`](https://mc-stan.org/loo/reference/loo.html) or
  [`waic`](https://mc-stan.org/loo/reference/waic.html), such as
  `k_threshold` or `cores`.

## Value

An object from the `loo` package: a `psis_loo` or `waic` object.
`loo_compare()` accepts several of them, which is how fits are compared;
the differences in ELPD are what matters, not the absolute value. For
`unit = "individual"`, the attribute `"units"` names the individual of
each pointwise term, so a high Pareto \\k\\ can be traced to an animal.

## Details

Relatives share breeding values, so an animal model does not factorize
over individuals and the density of one animal depends on all the
others. The default terms are those of leave-one-out for a
non-factorized normal model (Bürkner, Gabry and Vehtari 2021): for each
draw, the conditional density of an individual's observed traits given
every other individual, with that individual's own gaps integrated out.
When other individuals have gaps, the draw conditions on their sampled
values (`Y_mis`), which keeps the importance sampling exact. Individuals
with no observed trait carry no term. The cost per draw is two \\N
\times N\\ by \\N \times K\\ matrix products, batched over draws, plus
one decomposition of `A`.

Phenotypic models use one term per individual, integrating missing
traits out, from their own `log_lik`.

With `unit = "record"`, the terms are the animal model's `log_lik`: for
complete responses they sum to the joint log likelihood, but the
rotation mixes individuals, so a term is not one animal left out. That
criterion is refused for responses with gaps, whose `log_lik` terms use
imputed values.

Compare fits to the same response values, observations and missingness
pattern, with the same unit. The terms are densities of the standardized
response, so the absolute ELPD carries a constant that depends on the
data's scale; it cancels in comparisons between fits to the same
response. Individual-level terms assume the built-in likelihood and are
refused for a custom `model_file`. If response or model metadata are
unavailable (a bare draws matrix or an older fit), the model's `log_lik`
is used with a warning: the caller must verify that its terms are
suitable predictive units. Extracting raw draws and passing them
directly to [`loo::loo()`](https://mc-stan.org/loo/reference/loo.html)
or [`loo::waic()`](https://mc-stan.org/loo/reference/waic.html) bypasses
these checks.

**Fits to size and shape scores.** The scores of
[`size_shape`](https://diogro.github.io/stanimal/reference/size_shape.md)
are the same linear map of every individual's traits, \\s = M y + c\\.
Conditioning on the other individuals' scores is conditioning on their
traits, so a term of the score fit becomes the term of the traits by the
Jacobian \\\|\det M\|\\, a constant, together with the change from the
scores' standard deviations to the traits'. `transform` adds that
constant to every term, and `loo_compare()` can then set a fit to the
scores against a fit to the traits. The constant leaves the Pareto \\k\\
values and the standard error of a difference unchanged. It needs a
complete response (gaps are filled before the rotation, so the scores
are not a function of the observed data) and every shape axis (with
fewer, \\M\\ is not invertible). The shape axes are principal components
of the data, including the individual left out; that touches the prior
and not the likelihood, and with more than a few dozen individuals it is
negligible. The result carries `attr(, "scale") = "traits"`.

## References

Bürkner, P.-C., Gabry, J. and Vehtari, A. (2021). Efficient
leave-one-out cross-validation for Bayesian non-factorized normal and
Student-t models. *Computational Statistics* 36, 1243–1261.
[doi:10.1007/s00180-020-01045-4](https://doi.org/10.1007/s00180-020-01045-4)

## See also

[`posterior_predictive`](https://diogro.github.io/stanimal/reference/posterior_predictive.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit1 <- lmm_animal(Y, X, A)
fit2 <- lmm_animal(Y, X, A, genetic_structure = "factor", n_factors = 2)
loo::loo_compare(fit_loo(fit1), fit_loo(fit2))

# raw traits against size and shape scores
ss <- size_shape(Y)
fit3 <- lmm_animal(ss$scores, X, A)
loo::loo_compare(fit_loo(fit1), fit_loo(fit3, transform = ss))
} # }
```
