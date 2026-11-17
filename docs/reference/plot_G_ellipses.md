# Pairs plot of the data with posterior draws of G ellipses

For each pair of traits, the lower triangle shows the data – or their
residuals when `X` is given – with the `level` ellipses of a sample of
posterior G draws (thin lines), the posterior mean G (thick line) and
the covariance of those points (dashed), all centred on their mean. The
upper triangle reports the posterior mean genetic correlation and its
HPD interval.

## Usage

``` r
plot_G_ellipses(
  fit,
  Y,
  trait_names = NULL,
  n_show = 100,
  level = 0.95,
  prob = 0.95,
  X = NULL
)
```

## Arguments

- fit:

  A multivariate fit from `lmm_animal` (`cmd_stan = TRUE`) or a
  `posterior` draws object containing `G` and, when `X` is given, the
  fixed effects `beta`.

- Y:

  Matrix of the observed traits (the `Y` given to the fit). Missing
  values are allowed: the ellipses are built from the pairs of traits
  that were measured together, and a point is drawn only for the
  individuals that have both traits of its panel.

- trait_names:

  Optional trait names; default are the column names of `Y`.

- n_show:

  Number of posterior draws of G to draw, at least 1.

- level:

  Probability level of the ellipses, between 0 and 1.

- prob:

  Probability mass of the HPD interval of the correlations.

- X:

  Optional fixed-effects design matrix (the `X` given to the fit). When
  supplied, the traits are replaced by their residuals,
  `Y - X %*% t(beta)`, using the posterior mean of `beta`, and `P` is
  the covariance of those residuals. Pass it whenever the model has
  covariates or factors.

## Value

A ggplot object.

## Details

`G` and `E` describe the traits after the fixed effects, so with
`X = NULL` the points and the dashed ellipse include the sex, age or
group differences of the sample. Those effects can dwarf the residual
covariance, leaving a cloud much wider than the ellipses; supplying `X`
removes them and the scatter then has roughly the covariance
`P = G + E`. Residuals are formed with the posterior mean of `beta`, so
the points are those of a residual plot, not of an individual draw.
