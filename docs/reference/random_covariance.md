# Random covariance matrices

Draws random covariance matrices for simulations: a random correlation
matrix scaled by a vector of variances. Correlations come from the LKJ
distribution used as a prior by
[`lmm_animal`](https://diogro.github.io/stanimal/reference/lmm_animal.md),
from the partial-correlation construction of
[`evolqg::RandomMatrix()`](https://rdrr.io/pkg/evolqg/man/RandomMatrix.html),
or from a two-regime model that combines a positive size axis with a
sign-oriented random LKJ residual.

## Usage

``` r
random_covariance(
  n_traits,
  n_matrices = 1,
  variances = NULL,
  min_var = 1,
  max_var = 1,
  method = c("lkj", "partial", "factor"),
  lkj_prior = 2,
  integration = NULL,
  trait_names = NULL
)
```

## Arguments

- n_traits:

  Number of traits, a whole number of at least 1.

- n_matrices:

  Number of matrices to draw, a whole number of at least 1.

- variances:

  Optional vector of variances, one per trait, with finite positive
  values. The default draws them uniformly between `min_var` and
  `max_var`.

- min_var, max_var:

  Bounds of the uniform distribution of variances when `variances` is
  `NULL`: finite positive numbers with `min_var` not above `max_var`.

- method:

  Correlation model. `"lkj"` (default) draws a correlation matrix from
  the LKJ distribution with shape `lkj_prior`; `"partial"` draws partial
  correlations uniformly, reproducing
  [`evolqg::RandomMatrix()`](https://rdrr.io/pkg/evolqg/man/RandomMatrix.html);
  `"factor"` combines a positive size factor with a sign-oriented random
  LKJ residual, so PC1 is a size axis and draws differ at a fixed
  integration level.

- lkj_prior:

  LKJ shape for `method = "lkj"` and for the random residual of
  `method = "factor"`, a finite number above 0. The default `2` matches
  the prior of
  [`lmm_animal`](https://diogro.github.io/stanimal/reference/lmm_animal.md).

- integration:

  For `method = "factor"`, the target integration level, measured as the
  relative eigenvalue variance of `evolqg::CalcEigenVar(rel = TRUE)` and
  so between 0 and 1. A single number targets that level for every
  matrix; two numbers draw the target uniformly between them. The
  default `NULL` draws from 0.1 to 0.8, spanning weak to strong
  integration. Supplying it with another method is an error.

- trait_names:

  Optional trait names, one per trait; the default is `trait1`,
  `trait2`, ...

## Value

With `n_matrices = 1`, a named symmetric positive-definite matrix. With
more than one matrix, a list of such matrices.

## Details

The variances are independent of the correlations, and the correlation
matrix \\R\\ is scaled as \\R \odot \sqrt{v v^T}\\, so the diagonal is
\\v\\. With the default `min_var = max_var = 1` the result is a
correlation matrix, as in
[`evolqg::RandomMatrix()`](https://rdrr.io/pkg/evolqg/man/RandomMatrix.html);
using those defaults for both `G` and `E` gives a heritability of 0.5
for every trait. Use `variances`, or a wider `min_var`/`max_var`, to
vary the heritabilities.

The LKJ method is the default because it matches the prior family of the
models and stays well conditioned: correlations are shrunk toward zero
by `lkj_prior`, and a five-trait matrix typically has a condition number
in the tens. The partial-correlation method reproduces
[`evolqg::RandomMatrix()`](https://rdrr.io/pkg/evolqg/man/RandomMatrix.html)
and is the one to use when the simulation should span a realistic range
of integration: at five traits its first principal component can explain
up to about 90% of the variance and the mean squared correlation can
approach 0.8, while LKJ(2) stays well below that. The cost is
conditioning: the strongly integrated draws often have condition numbers
in the hundreds or more, whereas real size-dominated matrices have a
single positive axis and can reach the same integration with much
smaller condition numbers. Use the partial method to span the range, or
deliberately as a hard case for an estimator, not as a neutral default.

Unlike
[`evolqg::RandomMatrix()`](https://rdrr.io/pkg/evolqg/man/RandomMatrix.html),
this function also handles one and two traits with `method = "partial"`.

With `method = "factor"` the random part is \\R\_{res} = (1 - \delta)
I + \delta R\_{LKJ}\\, where \\R\_{LKJ}\\ is a random LKJ correlation
whose leading eigenvector has been oriented to be positive and \\\delta
= 0.3\\. Orienting keeps PC1 positive while leaving the remaining
components free to carry negative correlations. For targets below the
integration of \\R\_{res}\\, the matrix is \\(1 - \alpha) I + \alpha
R\_{res}\\, which shrinks the random structure toward independence
without flattening it into equicorrelation. For higher targets the
matrix is \\\rho 11^T + (1 - \rho) R\_{res}\\, which adds the
equal-loading size factor; at \\\rho = 1\\ it is all ones, whose
relative eigenvalue variance is 1. Both regimes have diagonal 1 and are
solved to match the requested integration. Because \\R\_{LKJ}\\ is
redrawn for every matrix, draws at a fixed integration differ in their
lower components. This is the method to use for size-dominated
morphometric simulations: it spans the full integration range while
staying far better conditioned than bending an LKJ spectrum directly (at
five traits the median condition number is about 4 at 0.1, 16 at 0.5 and
53 at 0.8). At two traits the target fixes the single correlation, so
draws at a fixed integration are identical up to the sign orientation.

## See also

[`animal_prior_draws`](https://diogro.github.io/stanimal/reference/animal_prior_draws.md),
[`size_shape`](https://diogro.github.io/stanimal/reference/size_shape.md)

## Examples

``` r
set.seed(1)
G <- random_covariance(3, min_var = 0.5, max_var = 1.5)
round(G, 2)
#>        trait1 trait2 trait3
#> trait1   0.77   0.51   0.19
#> trait2   0.51   0.87  -0.43
#> trait3   0.19  -0.43   1.07
round(cov2cor(G), 2)
#>        trait1 trait2 trait3
#> trait1   1.00   0.62   0.21
#> trait2   0.62   1.00  -0.44
#> trait3   0.21  -0.44   1.00
eigen(G, symmetric = TRUE, only.values = TRUE)$values
#> [1] 1.4772680 1.1197192 0.1134987

# Reproduce the evolqg construction
random_covariance(2, n_matrices = 3, method = "partial")
#> [[1]]
#>            trait1     trait2
#> trait1  1.0000000 -0.6468865
#> trait2 -0.6468865  1.0000000
#> 
#> [[2]]
#>           trait1    trait2
#> trait1 1.0000000 0.5396828
#> trait2 0.5396828 1.0000000
#> 
#> [[3]]
#>           trait1    trait2
#> trait1 1.0000000 0.9838122
#> trait2 0.9838122 1.0000000
#> 

# A strongly integrated size axis: PC1 dominates and is positive
S <- random_covariance(5, method = "factor", integration = 0.8)
R <- cov2cor(S)
round(R, 2)
#>        trait1 trait2 trait3 trait4 trait5
#> trait1   1.00   0.90   0.89   0.89   0.91
#> trait2   0.90   1.00   0.90   0.89   0.90
#> trait3   0.89   0.90   1.00   0.87   0.88
#> trait4   0.89   0.89   0.87   1.00   0.90
#> trait5   0.91   0.90   0.88   0.90   1.00
eigen(R, symmetric = TRUE)$vectors[, 1]
#> [1] -0.4493574 -0.4490820 -0.4438025 -0.4448958 -0.4488991
```
