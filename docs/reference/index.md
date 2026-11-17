# Package index

## Fitting models

- [`lmm_animal()`](https://diogro.github.io/stanimal/reference/lmm_animal.md)
  : Fit linear mixed animal model
- [`lm_phen_cov()`](https://diogro.github.io/stanimal/reference/lm_phen_cov.md)
  : Fit a Bayesian phenotypic covariance matrix
- [`lmm_animal_laplace()`](https://diogro.github.io/stanimal/reference/laplace.md)
  [`lm_phen_cov_laplace()`](https://diogro.github.io/stanimal/reference/laplace.md)
  : Laplace approximations of the covariance fits
- [`lmm_animal_gibbs()`](https://diogro.github.io/stanimal/reference/lmm_animal_gibbs.md)
  : Fit an animal model with a Gibbs sampler

## Building the inputs

- [`animal_design()`](https://diogro.github.io/stanimal/reference/animal_design.md)
  : Response and design matrices from a model formula
- [`makeA()`](https://diogro.github.io/stanimal/reference/makeA.md) :
  Additive relationship matrix from a pedigree
- [`genomic_relationship()`](https://diogro.github.io/stanimal/reference/genomic_relationship.md)
  : Genomic relationship matrix from genotype data
- [`size_shape()`](https://diogro.github.io/stanimal/reference/size_shape.md)
  : Size and shape scores for many correlated traits
- [`back_transform()`](https://diogro.github.io/stanimal/reference/back_transform.md)
  : Back-transform posterior draws to the original trait scale
- [`pc1_variance_share()`](https://diogro.github.io/stanimal/reference/pc1_variance_share.md)
  : Share of variance explained by the first principal component

## Choosing the number of factors

- [`suggest_n_factors()`](https://diogro.github.io/stanimal/reference/suggest_n_factors.md)
  : Number of factors suggested by parallel analysis
- [`scree_eigenvalues()`](https://diogro.github.io/stanimal/reference/scree_eigenvalues.md)
  : Eigenvalues of the residual trait correlation matrix
- [`plot(`*`<scree>`*`)`](https://diogro.github.io/stanimal/reference/plot.scree.md)
  : Scree plot

## Reading the posterior

- [`fit_draws()`](https://diogro.github.io/stanimal/reference/fit_draws.md)
  : Posterior draws of a fit, shaped by trait
- [`draws_mean()`](https://diogro.github.io/stanimal/reference/draws_mean.md)
  [`draws_quantile()`](https://diogro.github.io/stanimal/reference/draws_mean.md)
  : Posterior means or quantiles of a fit, as matrices
- [`rename_beta()`](https://diogro.github.io/stanimal/reference/rename_beta.md)
  : Posterior draws with the fixed effects named
- [`to_draws_list()`](https://diogro.github.io/stanimal/reference/to_draws_list.md)
  : Convert the back-transformed draws to posterior draws
- [`animal_parameter_draws()`](https://diogro.github.io/stanimal/reference/animal_parameter_draws.md)
  : Heritability and genetic correlation draws from a fit

## Post-processing a fit

- [`breeding_values()`](https://diogro.github.io/stanimal/reference/breeding_values.md)
  [`print(`*`<stanimal_breeding_values>`*`)`](https://diogro.github.io/stanimal/reference/breeding_values.md)
  : Estimated breeding values (BLUPs)
- [`imputed_data()`](https://diogro.github.io/stanimal/reference/imputed_data.md)
  [`print(`*`<stanimal_imputation>`*`)`](https://diogro.github.io/stanimal/reference/imputed_data.md)
  : Completed datasets from a fit with missing values
- [`fit_loo()`](https://diogro.github.io/stanimal/reference/fit_loo.md)
  [`loo(`*`<CmdStanMCMC>`*`)`](https://diogro.github.io/stanimal/reference/fit_loo.md)
  [`loo(`*`<stanfit>`*`)`](https://diogro.github.io/stanimal/reference/fit_loo.md)
  : Leave-one-out cross-validation for a fit
- [`posterior_predictive()`](https://diogro.github.io/stanimal/reference/posterior_predictive.md)
  [`print(`*`<stanimal_yrep>`*`)`](https://diogro.github.io/stanimal/reference/posterior_predictive.md)
  [`plot(`*`<stanimal_yrep>`*`)`](https://diogro.github.io/stanimal/reference/posterior_predictive.md)
  : Replicate the response from the posterior

## Plots

- [`plot_prior_posterior()`](https://diogro.github.io/stanimal/reference/plot_prior_posterior.md)
  [`plot_forest()`](https://diogro.github.io/stanimal/reference/plot_prior_posterior.md)
  : Posterior versus prior for heritabilities and genetic correlations
- [`plot_G_ellipses()`](https://diogro.github.io/stanimal/reference/plot_G_ellipses.md)
  : Pairs plot of the data with posterior draws of G ellipses
- [`plot_corr()`](https://diogro.github.io/stanimal/reference/plot_corr.md)
  : Heatmap of a correlation matrix

## Priors and simulation

- [`animal_prior_draws()`](https://diogro.github.io/stanimal/reference/animal_prior_draws.md)
  : Prior draws for heritabilities and genetic correlations
- [`random_covariance()`](https://diogro.github.io/stanimal/reference/random_covariance.md)
  : Random covariance matrices

## Data

- [`ratones_morpho`](https://diogro.github.io/stanimal/reference/ratones_morpho.md)
  : Morphological and life-history data for five mouse lines
- [`ratones_pedigree`](https://diogro.github.io/stanimal/reference/ratones_pedigree.md)
  : Pedigree of the mouse lines
- [`opossum`](https://diogro.github.io/stanimal/reference/opossum.md) :
  Opossum cranial measurements and a genomic relationship matrix

## Package

- [`stanimal-package`](https://diogro.github.io/stanimal/reference/stanimal-package.md)
  [`stanimal`](https://diogro.github.io/stanimal/reference/stanimal-package.md)
  : The 'stanimal' package.
- [`stan_cache_dir()`](https://diogro.github.io/stanimal/reference/stan_cache.md)
  [`stanimal_clear_cache()`](https://diogro.github.io/stanimal/reference/stan_cache.md)
  : Location of the compiled Stan model cache
