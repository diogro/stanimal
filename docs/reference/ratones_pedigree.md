# Pedigree of the mouse lines

The pedigree behind
[`ratones_morpho`](https://diogro.github.io/stanimal/reference/ratones_morpho.md):
all 525 individuals of the experiment, including the ancestors that were
never measured. Used to build the additive relationship matrix `A` for
[`lmm_animal`](https://diogro.github.io/stanimal/reference/lmm_animal.md).

## Usage

``` r
ratones_pedigree
```

## Format

A data frame with 525 rows and 3 columns, in the order expected by
[`animal_design`](https://diogro.github.io/stanimal/reference/animal_design.md):

- `id`:

  individual catalogue number;

- `sire`:

  catalogue number of the sire, `NA` for a founder or an unknown parent;

- `dam`:

  catalogue number of the dam, `NA` for a founder or an unknown parent.

## Source

Penna et al. (2017), archived at Dryad under CC0 1.0
([doi:10.5061/dryad.5gr8r](https://doi.org/10.5061/dryad.5gr8r) ) and
distributed by the MIT-licensed companion package
<https://github.com/diogro/ratones>.

## Details

The pedigree is complete and consistent, and the animal model will run
on it, but it does not make the mice a quantitative-genetic design. Read
the weak identification section of
[`ratones_morpho`](https://diogro.github.io/stanimal/reference/ratones_morpho.md)
before using it for anything beyond a worked example, and prefer
[`opossum`](https://diogro.github.io/stanimal/reference/opossum.md) and
its genomic relationship matrix `opossum_grm` when the genetic
covariance is the object of inference.

## References

Penna, A., Melo, D., Bernardi, S., Oyarzabal, M. I., & Marroig, G.
(2017). The evolution of phenotypic integration: How directional
selection reshapes covariation in mice. Evolution, 71(10), 2370-2380.
[doi:10.1111/evo.13304](https://doi.org/10.1111/evo.13304)

## See also

[`ratones_morpho`](https://diogro.github.io/stanimal/reference/ratones_morpho.md),
[`makeA`](https://diogro.github.io/stanimal/reference/makeA.md)

## Examples

``` r
data(ratones_pedigree)
dim(ratones_pedigree)
#> [1] 525   3
head(ratones_pedigree)
#>       id sire dam
#> 1 260482   NA  NA
#> 2 260386   NA  NA
#> 3 260460   NA  NA
#> 4 260466   NA  NA
#> 5 260474   NA  NA
#> 6 260396   NA  NA
sum(is.na(ratones_pedigree$sire))  # founders
#> [1] 84
```
