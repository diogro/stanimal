# Back-transform posterior draws to the original trait scale

Takes a fit of the size and shape scores made with
[`lmm_animal`](https://diogro.github.io/stanimal/reference/lmm_animal.md)
and returns draws of the genetic and residual covariance matrices, the
phenotypic covariance, the genetic correlations and the heritabilities
on the scale of the original traits.

## Usage

``` r
back_transform(fit, transform)
```

## Arguments

- fit:

  A fit returned by
  [`lmm_animal`](https://diogro.github.io/stanimal/reference/lmm_animal.md)
  (a `CmdStanMCMC` object), or a `posterior` draws object, or a matrix
  of draws, containing the variables `G` and `E`, or `P` for a
  phenotypic fit from
  [`lm_phen_cov`](https://diogro.github.io/stanimal/reference/lm_phen_cov.md).

- transform:

  A `size_shape` object from
  [`size_shape`](https://diogro.github.io/stanimal/reference/size_shape.md).

## Value

A list of class `size_shape_draws` holding one element per quantity, in
the same layout as
[`lmm_animal_gibbs`](https://diogro.github.io/stanimal/reference/lmm_animal_gibbs.md):
traits first and draws last, with every dimension named after the
traits. For an animal fit it holds `G`, `E`, `P`, `corrG` and `corrP` as
`traits` by `traits` by `draws` arrays, and `h2` as a `traits` by
`draws` matrix. For a phenotypic fit it holds `P` and `corrP` as arrays
and `sigma_P` as a matrix. `traits` holds the trait names.

## Details

For each draw the covariance matrices are mapped as \\M \to B M B^T\\,
which is exact because size and shape are a linear reparametrization of
the traits. Heritabilities and correlations are then recomputed from the
transformed matrices: unlike the covariance, they are not invariant
under a change of traits, so the h2 of a shape axis is not the h2 of any
original trait.

The fixed effects `beta` stay on the model scale, since a coefficient
for size or a shape axis does not correspond to a single original trait.

Draws are the last dimension, as in
[`lmm_animal_gibbs`](https://diogro.github.io/stanimal/reference/lmm_animal_gibbs.md),
so posterior summaries are direct: `draws_mean(post, "G")` is the
posterior mean of `G` on the trait scale, and
`draws_quantile(post, "corrG", probs = 0.05)` its 5% quantile.
[`fit_draws`](https://diogro.github.io/stanimal/reference/fit_draws.md)
also accepts this object and returns the same arrays.

## See also

[`size_shape`](https://diogro.github.io/stanimal/reference/size_shape.md),
[`rename_beta`](https://diogro.github.io/stanimal/reference/rename_beta.md),
[`to_draws_list`](https://diogro.github.io/stanimal/reference/to_draws_list.md)

## Examples

``` r
ss <- size_shape(cbind(a = rnorm(20), b = rnorm(20), c = rnorm(20)))
k <- length(ss$raw_traits)
names_mat <- function(var) sprintf(paste0(var, "[%d,%d]"),
                                   rep(seq_len(k), k), rep(seq_len(k), each = k))
# a real fit supplies G and E; here the same draw is repeated
G <- random_covariance(k, min_var = 0.5, max_var = 1.5)
E <- random_covariance(k, min_var = 0.5, max_var = 1.5)
m <- cbind(matrix(rep(as.vector(G), 4), 4, byrow = TRUE),
           matrix(rep(as.vector(E), 4), 4, byrow = TRUE))
colnames(m) <- c(names_mat("G"), names_mat("E"))
post <- back_transform(m, ss)
dim(post$G)                    # traits x traits x draws
#> [1] 3 3 4
draws_mean(post, "G")          # posterior mean G on the trait scale
#>           a         b         c
#> a 2.4148785 0.5599549 0.2250826
#> b 0.5599549 0.8963117 0.2172248
#> c 0.2250826 0.2172248 1.0309030
draws_mean(post, "h2")
#>         a         b         c 
#> 0.4483425 0.5086840 0.3443663 
```
