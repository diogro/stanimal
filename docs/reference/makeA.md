# Additive relationship matrix from a pedigree

Computes the numerator (additive) relationship matrix `A` from a
pedigree, for use as the `A` argument of
[`lmm_animal`](https://diogro.github.io/stanimal/reference/lmm_animal.md).

## Usage

``` r
makeA(pedigree, Gamma = NULL, groups = NULL)
```

## Arguments

- pedigree:

  Data frame whose first three columns are the individual, its sire and
  its dam, with `0`, `NA` or `""` for unknown parents. Individuals can
  be in any order, and the pedigree can contain ancestors that have no
  records.

- Gamma:

  Optional relationship matrix among the founders, or among founder
  groups when `groups` is given. Founders are the individuals with both
  parents unknown, including parents that were added because they were
  missing from the pedigree. The default `NULL` leaves founders
  unrelated, as in a standard pedigree. With `groups = NULL`, `Gamma`
  must have one row and column per founder, or founder IDs as row and
  column names; with `groups`, it must have one row and column per
  group, or the group labels as row and column names.

- groups:

  Optional vector assigning each founder to a row and column of `Gamma`,
  for a metafounder or genetic-group model. It can be named with founder
  IDs, or unnamed with one entry per founder in pedigree order
  (individuals with both parents unknown, then any parents added as
  founders). When `Gamma` has no dimnames, groups are used in
  `sort(unique(groups))` order.

## Value

A symmetric positive-definite matrix with the individuals' IDs on both
dimensions, in the order of `pedigree`.

## Details

The additive relationship is computed with the tabular method: with
\\s_i\\ and \\d_i\\ the parents of individual \\i\\, \$\$A\_{ii} = 1 +
A\_{s_i d_i} / 2, \qquad A\_{ij} = (A\_{s_i j} + A\_{d_i j}) / 2,\$\$
taking a term as zero when the parent is unknown. The first line
includes the inbreeding coefficient of \\i\\ through the relationship
between its parents, so pedigrees with inbreeding are handled.
Individuals are processed in an order where parents come before their
offspring, whatever order the pedigree is in, and a pedigree where an
individual is its own ancestor is an error.

A parent that is not one of the individuals of the pedigree is added to
`A` as a founder, with a warning, so a truncated pedigree or a mistyped
ID does not silently drop the relationship. Subset `A` by name to match
the rows of `Y`, or let
[`animal_design`](https://diogro.github.io/stanimal/reference/animal_design.md)
do it.

[`animal_design`](https://diogro.github.io/stanimal/reference/animal_design.md)
calls this function when it is given `pedigree` and no genotypes.

With `Gamma`, the founder block of `A` is initialized with the founder
relationships instead of the identity, and the same recursion propagates
them. This is the \\A^\Gamma\\ construction for related founders and
metafounders (Legarra et al. 2015): \\\Gamma\_{ii}\\ is the relationship
of founder \\i\\ with itself and \\\Gamma\_{ij}\\ its relationship with
founder \\j\\, so a founder no longer has to have \\A\_{ii} = 1\\. When
`groups` is given, founders in the same group share the diagonal of
`Gamma` and different groups use its off-diagonals. A genomic
relationship matrix among the founders, for example from
[`genomic_relationship`](https://diogro.github.io/stanimal/reference/genomic_relationship.md),
can be used as `Gamma`, but it estimates realized relationships and is
not the same as the base-population \\\Gamma\\; shrinkage, such as the
`blend` argument of `genomic_relationship`, is advisable.

Only founders with both parents unknown are given founder relationships.
An unknown parent of an individual that has one known parent is still
treated as an unrelated founder, as in the standard tabular method; to
model those unknown parents as part of the related base population, add
them explicitly to the pedigree and assign them to a group.

## References

Legarra, A., Christensen, O. F., Vitezica, Z. G., Aguilar, I., &
Misztal, I. (2015). Ancestral relationships using metafounders: finite
ancestral populations and across population relationships. Genetics,
200(2), 455-468.
[doi:10.1534/genetics.115.177014](https://doi.org/10.1534/genetics.115.177014)

## See also

[`genomic_relationship`](https://diogro.github.io/stanimal/reference/genomic_relationship.md),
[`animal_design`](https://diogro.github.io/stanimal/reference/animal_design.md)

## Examples

``` r
ped <- data.frame(id = 1:6, sire = c(0, 0, 1, 1, 4, 4),
                  dam = c(0, 0, 2, 2, 3, 5))
A <- makeA(ped)
round(A, 3)
#>     1   2     3     4    5     6
#> 1 1.0 0.0 0.500 0.500 0.50 0.500
#> 2 0.0 1.0 0.500 0.500 0.50 0.500
#> 3 0.5 0.5 1.000 0.500 0.75 0.625
#> 4 0.5 0.5 0.500 1.000 0.75 0.875
#> 5 0.5 0.5 0.750 0.750 1.25 1.000
#> 6 0.5 0.5 0.625 0.875 1.00 1.375
# full-sib mating: individual 6 is the offspring of the sibs 4 and 5
ped2 <- data.frame(id = 1:6, sire = c(NA, NA, NA, 1, 1, 4),
                   dam = c(NA, NA, NA, 2, 2, 5))
diag(makeA(ped2))   # 1.25 for the inbred individual
#>    1    2    3    4    5    6 
#> 1.00 1.00 1.00 1.00 1.00 1.25 

# related founders: founders 1 and 2 have relationship 0.5
Gamma <- matrix(c(1, 0.5, 0.5, 1), 2,
                dimnames = list(c("1", "2"), c("1", "2")))
ped3 <- data.frame(id = 1:4, sire = c(NA, NA, 1, 1),
                   dam = c(NA, NA, 2, 2))
A3 <- makeA(ped3, Gamma = Gamma)
round(A3, 3)
#>      1    2    3    4
#> 1 1.00 0.50 0.75 0.75
#> 2 0.50 1.00 0.75 0.75
#> 3 0.75 0.75 1.25 0.75
#> 4 0.75 0.75 0.75 1.25

# founder groups: 1 and 2 in group A, 3 in group B
Gamma_groups <- matrix(c(1, 0.2, 0.2, 1), 2,
                       dimnames = list(c("A", "B"), c("A", "B")))
ped4 <- data.frame(id = 1:5, sire = c(NA, NA, NA, 1, 1),
                   dam = c(NA, NA, NA, 2, 3))
round(makeA(ped4, Gamma = Gamma_groups, groups = c("A", "A", "B")), 3)
#>     1   2   3   4   5
#> 1 1.0 1.0 0.2 1.0 0.6
#> 2 1.0 1.0 0.2 1.0 0.6
#> 3 0.2 0.2 1.0 0.2 0.6
#> 4 1.0 1.0 0.2 1.5 0.6
#> 5 0.6 0.6 0.6 0.6 1.1
```
