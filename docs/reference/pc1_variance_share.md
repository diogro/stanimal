# Share of variance explained by the first principal component

Diagnostic for the size and shape reparametrization: when a single
component dominates, the traits are mostly one size gradient plus noise,
and fitting size and shape scores is a much better behaved model than
fitting all of the raw traits.

## Usage

``` r
pc1_variance_share(Y)
```

## Arguments

- Y:

  Numeric matrix of traits, individuals in rows. `NA` is allowed and is
  not discarded: the correlations are estimated from the available
  values, as described in
  [`size_shape`](https://diogro.github.io/stanimal/reference/size_shape.md).

## Value

The share of the total variance taken by the first principal component
of the correlation matrix of `Y`, between `1 / ncol(Y)` and 1 (1 for a
single trait).

## Details

The correlation matrix is used so that traits measured on different
scales do not dominate through their units alone. A share near one means
the traits are one gradient, which is what
[`size_shape`](https://diogro.github.io/stanimal/reference/size_shape.md)
is for; a share close to `1 / ncol(Y)` means the traits carry
independent information and no reparametrization is needed.

## See also

[`size_shape`](https://diogro.github.io/stanimal/reference/size_shape.md),
[`animal_design`](https://diogro.github.io/stanimal/reference/animal_design.md)

## Examples

``` r
set.seed(1)
G <- random_covariance(3, method = "factor", integration = 0.8)
Y <- matrix(rnorm(200 * 3), 200, 3) %*% chol(G)
pc1_variance_share(Y)                     # high: one dominant size axis
#> [1] 0.9200658
pc1_variance_share(matrix(rnorm(600), 200, 3))  # close to 1 / 3
#> [1] 0.3555601
```
