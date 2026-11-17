# Heritability and genetic correlation draws from a fit

Heritability and genetic correlation draws from a fit

## Usage

``` r
animal_parameter_draws(fit, trait_names = NULL)
```

## Arguments

- fit:

  A fit returned by `lmm_animal`, or a `posterior` draws object
  containing `h2` and, for several traits, `corrG`.

- trait_names:

  Optional trait names. Defaults to names recorded on the fit, or
  generic names when these are unavailable.

## Value

A long data frame with columns `type`, `parameter` and `value`.
