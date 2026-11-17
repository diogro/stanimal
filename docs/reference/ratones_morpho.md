# Morphological and life-history data for five mouse lines

Cranial measurements and individual information for mice from the
long-term artificial selection experiment on body weight analyzed by
Penna et al. (2017). Selection for body weight at 49 days of age was
carried out downwards in the `s` and `h` lines and upwards in the `sp`
and `hp` lines, alongside a randomly mated control line `t`; the
`s`/`sp` lines were inbred by limiting population size and the `h`/`hp`
lines by full-sib mating in the first generations.

## Usage

``` r
ratones_morpho
```

## Format

A data frame with 329 rows and 46 columns. Row names are the individual
IDs (`as.character(ID)`).

- `.id`:

  line label combining the selection direction and the line, e.g.
  `"control.t"`;

- `ID`:

  individual catalogue number;

- `SEX`:

  sex, `"F"` or `"M"`;

- `GER`:

  generation of the long-term experiment;

- `LIN`:

  experimental line: `t` (control), `s` and `h` (downward selection),
  `sp` and `hp` (upward selection);

- `MADRE`:

  catalogue number of the dam;

- `PADRE`:

  catalogue number of the sire;

- `AGE`:

  age at sacrifice, in days;

- `P49`:

  body weight at 49 days of age;

- `line`:

  the line as a factor with levels `t`, `h'}, \code{s'`, `h` and `s`;

- `selection`:

  direction of selection: `"control"`, `"downwards"` or `"upwards"`;

- `IS_PM`, `IS_NSL`, `IS_PNS`, `PM_ZS`, `PM_ZI`, `PM_MT`, `NSL_NA`,
  `NSL_ZS`, `NSL_ZI`, `NA_BR`, `NA_PNS`, `BR_PT`, `BR_APET`, `PT_APET`,
  `PT_BA`, `PT_EAM`, `PT_ZYGO`, `PT_TSP`, `ZS_ZI`, `ZI_MT`, `ZI_ZYGO`,
  `ZI_TSP`, `MT_PNS`, `PNS_APET`, `APET_BA`, `APET_TS`, `BA_EAM`,
  `EAM_ZYGO`, `ZYGO_TSP`, `LD_AS`, `BR_LD`, `OPI_LD`, `PT_AS`, `JP_AS`,
  `BA_OPI`:

  35 cranial Euclidean distances, named after the landmarks they span.
  They are on the measured scale;
  [`log()`](https://rdrr.io/r/base/Log.html) is the usual transformation
  before fitting.

## Source

Penna et al. (2017), archived at Dryad under CC0 1.0
([doi:10.5061/dryad.5gr8r](https://doi.org/10.5061/dryad.5gr8r) ) and
distributed by the MIT-licensed companion package
<https://github.com/diogro/ratones>.

## Details

Every individual carries a full set of 35 cranial Euclidean distances,
measured as the mean of two replicas, plus the catalogue numbers of its
parents. The distances are moderately integrated: the first principal
component of their correlation matrix explains about 40% of the variance
and the mean pairwise correlation is 0.35, so the traits carry partly
independent information and a handful of them is enough for a
well-identified model.

## Weak identification

This is a long-term artificial-selection experiment that happens to
carry a pedigree, not a design built to estimate a genetic covariance
matrix. The animal model fits, but 329 individuals from 58 sires
identify a 35-trait \\\mathbf{G}\\ only weakly: the split between the
genetic and the environmental covariance is left largely to the priors,
so a heritability read from
[`lmm_animal`](https://diogro.github.io/stanimal/reference/lmm_animal.md)
on these traits depends on the prior as much as on the data, and can
move by a factor of eight between defensible choices. Use it as the
small worked example it is, and use
[`opossum`](https://diogro.github.io/stanimal/reference/opossum.md),
which carries marker genotypes and a design that does identify
\\\mathbf{G}\\, for anything about genetic covariance.

## References

Penna, A., Melo, D., Bernardi, S., Oyarzabal, M. I., & Marroig, G.
(2017). The evolution of phenotypic integration: How directional
selection reshapes covariation in mice. Evolution, 71(10), 2370-2380.
[doi:10.1111/evo.13304](https://doi.org/10.1111/evo.13304)

## See also

[`ratones_pedigree`](https://diogro.github.io/stanimal/reference/ratones_pedigree.md),
[`animal_design`](https://diogro.github.io/stanimal/reference/animal_design.md),
[`lm_phen_cov`](https://diogro.github.io/stanimal/reference/lm_phen_cov.md)

## Examples

``` r
data(ratones_morpho)
dim(ratones_morpho)
#> [1] 329  46
table(ratones_morpho$selection)
#> 
#>   control downwards   upwards 
#>        67       121       141 
traits <- names(ratones_morpho)[12:46]
pc1_variance_share(log(as.matrix(ratones_morpho[, traits])))
#> [1] 0.4063927
```
