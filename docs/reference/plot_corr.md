# Heatmap of a correlation matrix

Heatmap of a correlation matrix

## Usage

``` r
plot_corr(M, title = "")
```

## Arguments

- M:

  Correlation matrix, with optional row and column names.

- title:

  Plot title.

## Value

A ggplot object with a diverging red-blue scale fixed to \[-1, 1\].

## Examples

``` r
plot_corr(stats::cor(mtcars[, 1:5]), "mtcars")
```
