# Scree plot

Scree plot

## Usage

``` r
# S3 method for class 'scree'
plot(x, reference = NULL, title = "Scree plot", ...)
```

## Arguments

- x:

  A `scree` object from
  [`scree_eigenvalues`](https://diogro.github.io/stanimal/reference/scree_eigenvalues.md).

- reference:

  Optional eigenvalue cutoffs, either the object returned by
  [`suggest_n_factors`](https://diogro.github.io/stanimal/reference/suggest_n_factors.md)
  or its `reference` attribute. Components above the dashed line are the
  ones parallel analysis keeps.

- title:

  Plot title.

- ...:

  Unused, for consistency with `plot`.

## Value

A ggplot object.

## See also

[`scree_eigenvalues`](https://diogro.github.io/stanimal/reference/scree_eigenvalues.md),
[`suggest_n_factors`](https://diogro.github.io/stanimal/reference/suggest_n_factors.md)
