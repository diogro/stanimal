# Response and design matrices from a model formula

Builds the `Y` and `X` arguments of
[`lmm_animal`](https://diogro.github.io/stanimal/reference/lmm_animal.md),
[`lm_phen_cov`](https://diogro.github.io/stanimal/reference/lm_phen_cov.md)
and
[`lmm_animal_gibbs`](https://diogro.github.io/stanimal/reference/lmm_animal_gibbs.md)
from a formula and a data frame, so a model written the usual way
(`y ~ sex + age`, or `cbind(x, y, z) ~ sex + age`) can be fitted without
building the matrices by hand.

## Usage

``` r
animal_design(
  formula,
  data,
  na.action = stats::na.omit,
  pedigree = NULL,
  genotypes = NULL,
  Gamma = NULL,
  groups = NULL,
  ...
)
```

## Arguments

- formula:

  Model formula. The left-hand side is a single response or
  [`cbind()`](https://rdrr.io/r/base/cbind.html) of several traits; the
  right-hand side lists the fixed effects.

- data:

  Data frame holding the variables in `formula`. When `pedigree` or
  `genotypes` are used, its row names must be the individual IDs
  (`rownames(data) <- data$id`); an error says so when they were left to
  R.

- na.action:

  How to handle rows with missing *predictors*, passed to
  [`model.frame`](https://rdrr.io/r/stats/model.frame.html); the default
  drops them, since a row with an unknown fixed effect cannot be placed
  in the design. Missing values in the response are kept: a trait that
  was not measured on an individual is information the model can use,
  and
  [`lmm_animal()`](https://diogro.github.io/stanimal/reference/lmm_animal.md)
  integrates it out.

- pedigree:

  Optional pedigree, used to compute the additive relationship matrix
  `A` with
  [`makeA`](https://diogro.github.io/stanimal/reference/makeA.md). A
  data frame whose first three columns are the individual, its sire and
  its dam, with `0`, `NA` or `""` for unknown parents. The pedigree may
  contain ancestors that are not in `data`.

- genotypes:

  Optional genotype matrix, used to compute a genomic relationship
  matrix (see
  [`genomic_relationship`](https://diogro.github.io/stanimal/reference/genomic_relationship.md)):
  individuals in rows, with their IDs as row names, and biallelic
  markers in columns, coded 0, 1, 2.

- Gamma:

  Optional relationship matrix among pedigree founders, or among founder
  groups when `groups` is given; passed to
  [`makeA`](https://diogro.github.io/stanimal/reference/makeA.md). Only
  used with `pedigree`.

- groups:

  Optional founder-to-group assignment passed to
  [`makeA`](https://diogro.github.io/stanimal/reference/makeA.md); only
  used with `pedigree` and `Gamma`.

- ...:

  Further arguments passed to
  [`genomic_relationship`](https://diogro.github.io/stanimal/reference/genomic_relationship.md),
  such as `p`, `maf` or `blend`.

## Value

A list with

- Y:

  response matrix, one column per trait, named after the response(s) in
  the formula;

- X:

  fixed-effects design matrix, with columns named after the terms in the
  formula;

- A:

  relationship matrix, present only when `pedigree` or `genotypes` is
  given, restricted to the individuals in `Y` and in the same row and
  column order.

Rows are in the order of `data`. `attr(, "individuals")` holds the row
names of `data` for the rows kept, so naming the rows with the
individual IDs lets a relationship matrix `A` be matched to `Y`:
[`lmm_animal`](https://diogro.github.io/stanimal/reference/lmm_animal.md)
puts a named `A` in the row order of `Y` and errors if the names
disagree.

## Details

Random-effect terms such as `(1 | animal)` are ignored, so formulas
written for `lme4` can be reused unchanged: relatedness is supplied
through `A`, which is passed to the model directly.

The result satisfies what the model functions expect: `Y` is numeric and
`X` has full column rank. Factor predictors are dummy coded and the
intercept is kept unless the formula removes it, so the columns of `X`
match the `beta` draws reported by a fit. Removing the intercept is only
allowed when the remaining columns still sum to a constant, as the 0/1
indicators of a factor do, because the models centre the responses and
report their mean back through the design; a formula such as `y ~ 0 + x`
is an error.

When `pedigree` or `genotypes` is given, `A` is computed for all
individuals in the pedigree or genotype matrix and then subset to the
rows of `Y`, so ancestors without records and individuals with missing
phenotypes are handled: every row of `Y` must appear in the pedigree or
in the row names of the genotype matrix. With `pedigree`, `Gamma` and
`groups` allow related founders or metafounder groups; see
[`makeA`](https://diogro.github.io/stanimal/reference/makeA.md).

If the response has three or more traits and its first principal
component explains at least 80% of the variance, a warning suggests
fitting size and shape scores instead (see
[`size_shape`](https://diogro.github.io/stanimal/reference/size_shape.md));
[`pc1_variance_share`](https://diogro.github.io/stanimal/reference/pc1_variance_share.md)
reports that share. If the diagnostic cannot be computed, a warning
explains why; the design is still returned.

## See also

[`lmm_animal`](https://diogro.github.io/stanimal/reference/lmm_animal.md),
[`makeA`](https://diogro.github.io/stanimal/reference/makeA.md),
[`genomic_relationship`](https://diogro.github.io/stanimal/reference/genomic_relationship.md),
[`rename_beta`](https://diogro.github.io/stanimal/reference/rename_beta.md)

## Examples

``` r
dat <- data.frame(trait = rnorm(6), x = rnorm(6),
                  sex = rep(c(0, 1), each = 3), id = 1:6)
rownames(dat) <- dat$id
design <- animal_design(trait ~ sex + x, dat)
design$Y
#>          trait
#> 1 -1.400043517
#> 2  0.255317055
#> 3 -2.437263611
#> 4 -0.005571287
#> 5  0.621552721
#> 6  1.148411606
design$X
#>   (Intercept) sex          x
#> 1           1   0 -1.8218177
#> 2           1   0 -0.2473253
#> 3           1   0 -0.2441996
#> 4           1   1 -0.2827054
#> 5           1   1 -0.5536994
#> 6           1   1  0.6289820
#> attr(,"assign")
#> [1] 0 1 2

# formulas written for lme4 work as they are; the random term is ignored
animal_design(trait ~ sex + (1 | id), dat)
#> $Y
#>          trait
#> 1 -1.400043517
#> 2  0.255317055
#> 3 -2.437263611
#> 4 -0.005571287
#> 5  0.621552721
#> 6  1.148411606
#> 
#> $X
#>   (Intercept) sex
#> 1           1   0
#> 2           1   0
#> 3           1   0
#> 4           1   1
#> 5           1   1
#> 6           1   1
#> attr(,"assign")
#> [1] 0 1
#> 
#> attr(,"individuals")
#> [1] "1" "2" "3" "4" "5" "6"

# relatedness from a pedigree (first three columns: id, sire, dam)
ped <- data.frame(id = 1:6, sire = c(0, 0, 1, 1, 4, 4),
                  dam = c(0, 0, 2, 2, 3, 5))
design <- animal_design(trait ~ sex, dat, pedigree = ped)
round(design$A, 2)
#>     1   2    3    4    5    6
#> 1 1.0 0.0 0.50 0.50 0.50 0.50
#> 2 0.0 1.0 0.50 0.50 0.50 0.50
#> 3 0.5 0.5 1.00 0.50 0.75 0.62
#> 4 0.5 0.5 0.50 1.00 0.75 0.88
#> 5 0.5 0.5 0.75 0.75 1.25 1.00
#> 6 0.5 0.5 0.62 0.88 1.00 1.38

# related founders: pass their relationship matrix as Gamma
Gamma <- matrix(c(1, 0.5, 0.5, 1), 2,
                dimnames = list(c("1", "2"), c("1", "2")))
design <- animal_design(trait ~ sex, dat, pedigree = ped, Gamma = Gamma)
round(design$A, 2)
#>      1    2    3    4    5    6
#> 1 1.00 0.50 0.75 0.75 0.75 0.75
#> 2 0.50 1.00 0.75 0.75 0.75 0.75
#> 3 0.75 0.75 1.25 0.75 1.00 0.88
#> 4 0.75 0.75 0.75 1.25 1.00 1.12
#> 5 0.75 0.75 1.00 1.00 1.38 1.19
#> 6 0.75 0.75 0.88 1.12 1.19 1.50

# relatedness from genotypes (individuals in rows, markers in columns)
M <- matrix(rbinom(6 * 40, 2, 0.3), 6,
            dimnames = list(as.character(1:6), NULL))
design <- animal_design(trait ~ sex, dat, genotypes = M, blend = 0.01)
dim(design$A)
#> [1] 6 6
```
