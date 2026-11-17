# Opossum cranial measurements and a genomic relationship matrix

Cranial Euclidean distances and individual information for opossums from
the colony of the South Texas Diabetes and Obesity Institute, analysed
by Porto et al. (2016). The traits are the same 35 landmarks as
[`ratones_morpho`](https://diogro.github.io/stanimal/reference/ratones_morpho.md),
so the two datasets can be compared directly.

## Usage

``` r
opossum
```

## Format

Two objects.

- `opossum`:

  a data frame with 553 rows and 38 columns: `id`, `sex`, `age` in days,
  and the 35 cranial distances as measured, named after the landmarks
  they span, with underscores as in
  [`ratones_morpho`](https://diogro.github.io/stanimal/reference/ratones_morpho.md).
  Row names are the individual IDs. 22 individuals are missing every
  trait.
  [`lmm_animal`](https://diogro.github.io/stanimal/reference/lmm_animal.md)
  accepts fully missing response rows, but they must be removed before
  fitting
  [`lm_phen_cov`](https://diogro.github.io/stanimal/reference/lm_phen_cov.md).
  505 individuals have no missing value in any trait or covariate.

- `opossum_grm`:

  a 553 by 553 genomic relationship matrix computed from 3,696 markers
  as \\ZZ^\top / (2\sum_j p_j(1-p_j))\\ (VanRaden 2008), \\Z\\ being the
  allele-count matrix centred by twice the sample allele frequencies
  \\p_j\\. It equals
  [`genomic_relationship()`](https://diogro.github.io/stanimal/reference/genomic_relationship.md)
  applied to the markers, and it is on the scale of a pedigree
  relationship matrix, so `G` fitted on it is the additive genetic
  covariance of the population. Rows and columns are named with the
  individual IDs of `opossum`, so the rows of the subset being analysed
  can be passed as the `A` argument of
  [`lmm_animal`](https://diogro.github.io/stanimal/reference/lmm_animal.md)
  directly, as in the example. The full matrix is only positive
  semidefinite, so pass `blend_A = TRUE` if it is given whole. The
  diagonal runs from 0.64 to 1.70 with a mean of 1.05, and the
  off-diagonal entries average zero, since the allele frequencies come
  from the same animals. Versions up to 0.8.1 shipped the same matrix
  divided by the number of markers instead, 2.86 times smaller, which
  inflated `G` and the heritabilities fitted on it.

## Source

Porto, A., Schmelter, R., VandeBerg, J. L., Marroig, G., & Cheverud, J.
M. (2016). Evolution of the genotype-to-phenotype map and the cost of
pleiotropy in mammals. Genetics, 204(4), 1601-1612.
[doi:10.1534/genetics.116.189431](https://doi.org/10.1534/genetics.116.189431)
. Data deposited on figshare as “Data - Genetics - Porto et al 2016”,
[doi:10.6084/m9.figshare.4055961](https://doi.org/10.6084/m9.figshare.4055961)
.

The coded marker files in that archive carry no individual identifiers,
so the animals are matched to the cranial measurements by row order. The
two sources agree in length at 553, and the identities are attached here
as the row and column names of `opossum_grm`.

## Details

Unlike the mice, these animals carry marker genotypes, so the
relationship matrix is genomic rather than genealogical. That matters
for these data: the colony pedigree shipped with the original archive
gives every pair a kinship coefficient between 0.69 and 0.93, which
leaves additive variance indistinguishable from a shared colony effect.
The genomic matrix separates relatives from non-relatives and makes the
genetic covariance identifiable.

The distances are stored as measured. Porto et al. (2016) analysed them
untransformed, and the phenomenon is confirmed by their published
phenotypic covariance, whose diagonal correlates 0.834 with the
variances of the raw distances against 0.046 with those of their
logarithms.

## References

Porto, A., Schmelter, R., VandeBerg, J. L., Marroig, G., & Cheverud, J.
M. (2016). Evolution of the genotype-to-phenotype map and the cost of
pleiotropy in mammals. Genetics, 204(4), 1601-1612.
[doi:10.1534/genetics.116.189431](https://doi.org/10.1534/genetics.116.189431)

## See also

[`ratones_morpho`](https://diogro.github.io/stanimal/reference/ratones_morpho.md),
[`lmm_animal`](https://diogro.github.io/stanimal/reference/lmm_animal.md)

## Examples

``` r
data(opossum)
dim(opossum)
#> [1] 553  38
dim(opossum_grm)
#> [1] 553 553
traits <- setdiff(names(opossum), c("id", "sex", "age"))
if (FALSE) { # \dontrun{
# Phenotypic covariance of four traits, no pedigree needed
lm_phen_cov(as.matrix(opossum[, traits[1:4]]),
            model.matrix(~ sex + age, opossum))

# Genomic animal model on the same four traits
design <- animal_design(
  stats::reformulate(c("sex", "age"),
                     response = sprintf("cbind(%s)", paste(traits[1:4],
                                                           collapse = ", "))),
  data = opossum, pedigree = NULL)
lmm_animal(design$Y, design$X, opossum_grm[rownames(design$Y),
                                           rownames(design$Y)])
} # }
```
