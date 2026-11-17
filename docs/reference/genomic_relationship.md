# Genomic relationship matrix from genotype data

Computes the VanRaden (2008) genomic relationship matrix (GRM) from
biallelic marker genotypes, for use as the `A` argument of
[`lmm_animal`](https://diogro.github.io/stanimal/reference/lmm_animal.md).

## Usage

``` r
genomic_relationship(M, p = NULL, maf = 0.01, max_missing = 0.1, blend = 0)
```

## Arguments

- M:

  Numeric genotype matrix with individuals in rows and markers in
  columns, coded as the number of copies (0, 1, 2) of the counted
  allele. Missing genotypes are allowed (`NA`); non-integer values such
  as imputed dosages are also accepted.

- p:

  Optional vector of allele frequencies of the counted allele, one per
  column of `M`. If `NULL` (default), they are estimated from `M` as
  `colMeans(M, na.rm = TRUE) / 2`. Supplying base-population frequencies
  is recommended when they are known.

- maf:

  Markers with minor allele frequency below this value are dropped.
  Monomorphic markers (and, when `p` is estimated, markers with a single
  observed genotype class) are always dropped.

- max_missing:

  Markers with a proportion of missing genotypes above this value are
  dropped. Fully missing markers are always dropped, even when
  `max_missing = 1` and `p` is supplied.

- blend:

  Weight `w` in `(1 - w) * G + w * I`. Default 0 returns the raw GRM; a
  small value (e.g. 0.01) makes a singular GRM, such as one with fewer
  markers than individuals or duplicated genotypes, positive definite.

## Value

A symmetric positive-semidefinite matrix with the individuals' row names
of `M` on both dimensions. The number of markers used and their scaling
constant are stored in the attributes `n_markers` and `scale`.

## Details

The GRM is `G = Z %*% t(Z) / (2 * sum(p * (1 - p)))`, where
`Z = M - 2p`. Missing genotypes are replaced by their expected value
`2p`, so they contribute zero to the numerator. Fully missing markers
are excluded from both numerator and denominator; if other markers
remain, a warning is issued when they would otherwise pass the filters.
With frequencies estimated from the data, the diagonal averages about
one plus the inbreeding coefficient, and the off-diagonals are relative
to the sample average.

## References

VanRaden, P.M. (2008). Efficient methods to compute genomic predictions.
Journal of Dairy Science 91: 4414-4423.

## Examples

``` r
set.seed(1)
M <- matrix(rbinom(50 * 500, 2, 0.3), nrow = 50)
G <- genomic_relationship(M, blend = 0.01)
dim(G)
#> [1] 50 50
```
