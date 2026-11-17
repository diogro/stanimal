# Size and shape scores for many correlated traits

Rewrites a set of correlated traits as log size plus orthogonal shape
axes, the reparametrization used to fit many allometric traits: raw
traits are usually dominated by a single size gradient, so their genetic
covariance is a long, thin ridge that is slow to sample. Size and shape
span the same space, so nothing is lost, but the sampler works on
well-separated axes.

## Usage

``` r
size_shape(Y, n_shape = NULL)
```

## Arguments

- Y:

  Numeric matrix of traits, individuals in rows and traits in columns.
  Log-transformed measurements are the usual input. Missing values
  (`NA`) are allowed and are not discarded: see Details.

- n_shape:

  Number of shape axes to keep, between 1 and `ncol(Y) - 1`, and no more
  than the numerical rank of the centered shape data. The default
  requests all `ncol(Y) - 1` axes and errors if the data cannot identify
  them; choose fewer explicitly in that case. Keeping all axes makes the
  transformation a bijection; fewer axes give a lower-rank (singular)
  covariance in the original trait space.

## Value

A list of class `size_shape` with

- scores:

  matrix with the size and shape scores, one column each, named `size`,
  `shape1`, `shape2`, ...;

- loadings:

  `ncol(Y)` by `1 + n_shape` matrix `B` such that, up to a constant
  absorbed by the model intercepts, `Y = scores %*% t(B)`;

- traits:

  names of the model traits (`size`, `shape1`, ...);

- raw_traits:

  names of the original traits;

- centre:

  column means removed from the shape axes;

- gaps:

  number of entries that were not measured, and so were completed before
  the transformation (`0` for a complete response);

- scores_sd:

  standard deviation of each completed score, on the same layout as
  `scores`, or `NULL` when nothing was completed. Observed individuals
  have `0`.

## Details

Size is the mean of the traits, `rowMeans(Y)`, i.e. the log geometric
mean for log measurements. Shape is `Y - size`, which sums to zero
within each individual and therefore lives in an `ncol(Y) - 1`
dimensional subspace; its principal components are the shape axes. Size
is often close to the leading axis of variation of a morphometric data
set, but it need not be: it is the unweighted mean of the traits, while
the principal components are shaped by their covariance structure. The
shape axes are ordered by the residual variance they explain. Their
numerical rank is at most `min(nrow(Y) - 1, ncol(Y) - 1)` and can be
smaller if the shape data are collinear. It is computed after completing
missing values, using a relative singular-value tolerance of
`max(dim(Y)) * .Machine$double.eps`.

The transformation is linear, so genetic and residual covariance
matrices transform as \\B M B^T\\; use
[`back_transform`](https://diogro.github.io/stanimal/reference/back_transform.md)
to do this for the draws of a fit. Signs of the shape axes are
arbitrary.

The prior is not invariant under this change of traits, and that is
usually an asset. The LKJ prior that
[`lmm_animal`](https://diogro.github.io/stanimal/reference/lmm_animal.md)
puts on the correlation matrix acts on correlations between size and
shape instead of correlations between raw traits. The shape axes are
principal components of `Y - size`, so they are uncorrelated in the
sample and their genetic correlations are expected to be small; a larger
`lkj_prior` than for the same traits fitted raw (stronger shrinkage
towards zero correlation) is therefore reasonable, and it can also make
the model cheaper to sample.

**Missing values.** Size is a mean of all the traits and shape is a
difference between them, so a row with a gap has neither. The gap is
therefore filled first: the mean and covariance of the traits are
estimated by maximum likelihood from the values that were measured (an
EM fit of a multivariate normal, which needs no other information), and
a missing entry is replaced by its conditional mean given the traits
that individual does have. The rest of the transformation then runs
unchanged, which is exact for the observed values and a conditional
estimate for the completed ones.

Because those completed scores are estimates, `scores_sd` reports the
standard deviation each one carries, and it is worth reading before
comparing a mouse that was missing a trait with one that was not: an
entry completed from a single measured trait is mostly the model's
prediction. A row with nothing observed carries no information and is an
error, and two traits never measured on the same individual make their
covariance unidentifiable, which is also an error. The EM fit starts
from a positive-definite diagonal covariance and errors if it fails to
converge or produces a singular covariance. Check the missingness
pattern, collinearity and sample size in that case.

The completion is a point estimate, not a multiple imputation: the
scores are then treated as data by the model they are passed to, and the
uncertainty in `scores_sd` does not travel into that fit. It is the
transformation that is estimated once, exactly as it is for a complete
response.

## See also

[`back_transform`](https://diogro.github.io/stanimal/reference/back_transform.md),
[`animal_design`](https://diogro.github.io/stanimal/reference/animal_design.md),
[`imputed_data`](https://diogro.github.io/stanimal/reference/imputed_data.md)

## Examples

``` r
set.seed(1)
G <- random_covariance(3, method = "factor", integration = 0.8)
Y <- matrix(rnorm(30 * 3), 30, 3) %*% chol(G)
colnames(Y) <- c("a", "b", "c")
ss <- size_shape(Y)
head(ss$scores)
#>            size       shape1      shape2
#> [1,]  1.4713091 -0.174043125  0.09388516
#> [2,] -0.1054167 -0.527625904  0.23599427
#> [3,] -0.8306718 -0.001639013  0.13899487
#> [4,]  0.1249348 -0.625704797 -0.20420978
#> [5,]  0.8642806  0.270603383  0.34168101
#> [6,]  0.8045399  0.168324928 -0.19979879
ss$loadings
#>   size        PC1        PC2
#> a    1 -0.7464609  0.3308517
#> b    1  0.0867045 -0.8118799
#> c    1  0.6597564  0.4810283

# a gap in one trait is completed from the other two
Y[2, "b"] <- NA
size_shape(Y)$scores[2, ]
#>         size       shape1       shape2 
#>  0.017864835 -0.499582949 -0.000827619 
```
