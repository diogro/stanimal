# Convert the back-transformed draws to posterior draws

Reshapes the arrays returned by
[`back_transform`](https://diogro.github.io/stanimal/reference/back_transform.md)
into
[`draws_matrix`](https://mc-stan.org/posterior/reference/draws_matrix.html)
objects, one per quantity, so that the posterior and bayesplot tools can
be used on them.

## Usage

``` r
to_draws_list(x)
```

## Arguments

- x:

  A `size_shape_draws` object from
  [`back_transform`](https://diogro.github.io/stanimal/reference/back_transform.md).

## Value

A named list of `draws_matrix` objects, one per quantity of `x`, with
the variables named as the models report them: `G[i,j]`, `E[i,j]` and
`P[i,j]` for every element of the covariance matrices, `corrG[i,j]` and
`corrP[i,j]` for the correlations above the diagonal only, and `h2[i]`
or `sigma_P[i]` for the vectors. The trait names are kept in
`attr(, "traits")`.

## Details

Correlations are reported once per pair, for `i < j`, because the two
halves of the matrix are identical and a draws object has no notion of
symmetry.

To get everything back in a single draws object, bind the elements with
`do.call(posterior::bind_draws, unname(draws))`. The names of the list
have to be dropped, because `bind_draws()` takes the draws objects as
unnamed arguments.

## See also

[`back_transform`](https://diogro.github.io/stanimal/reference/back_transform.md)

## Examples

``` r
ss <- size_shape(cbind(a = rnorm(20), b = rnorm(20), c = rnorm(20)))
k <- length(ss$raw_traits)
G <- random_covariance(k, min_var = 0.5, max_var = 1.5)
E <- random_covariance(k, min_var = 0.5, max_var = 1.5)
m <- cbind(matrix(rep(as.vector(G), 3), 3, byrow = TRUE),
           matrix(rep(as.vector(E), 3), 3, byrow = TRUE))
colnames(m) <- c(sprintf("G[%d,%d]", rep(seq_len(k), k), rep(seq_len(k), each = k)),
                 sprintf("E[%d,%d]", rep(seq_len(k), k), rep(seq_len(k), each = k)))
draws <- to_draws_list(back_transform(m, ss))
names(draws)
#> [1] "G"     "E"     "P"     "corrG" "corrP" "h2"   
posterior::summarise_draws(draws$corrG)
#> # A tibble: 3 × 10
#>   variable    mean median    sd   mad    q5   q95  rhat ess_bulk ess_tail
#>   <chr>      <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
#> 1 corrG[1,2] 0.496  0.496     0     0 0.496 0.496    NA       NA       NA
#> 2 corrG[1,3] 0.199  0.199     0     0 0.199 0.199    NA       NA       NA
#> 3 corrG[2,3] 0.243  0.243     0     0 0.243 0.243    NA       NA       NA
combined <- do.call(posterior::bind_draws, unname(draws))
posterior::variables(combined)[1:4]
#> [1] "G[1,1]" "G[2,1]" "G[3,1]" "G[1,2]"
```
