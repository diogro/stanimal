skip_if_no_cmdstan()

phenotypic_factor_truth <- function() {
  lambda <- matrix(c(0.8, 0.0,
                     0.7, 0.4,
                     0.5, 0.0), 3, 2, byrow = TRUE)
  psi <- c(0.36, 0.40, 0.45)
  lambda %*% t(lambda) + diag(psi)
}

test_that("factor P model samples and exposes factor summaries", {
  P <- phenotypic_factor_truth()
  sim <- simulate_animal(n_fam = 1, size = 300, E = P, seed = 51,
                         genetic = FALSE)
  fit <- suppressWarnings(lm_phen_cov(
    sim$Y, sim$X, covariance_structure = "factor",
    n_factors = 2, hs_global_scale = 2, hs_slab_scale = 2,
    column_decay = 0.9, iter = 1000, warmup = 500,
    chains = 2, cores = 2, seed = 1, refresh = 0, show_messages = FALSE
  ))
  expect_s3_class(fit, "CmdStanMCMC")
  expect_clean_fit(fit)
  expect_lt(max(abs(post_matrix(fit, "P", 3) - P)), 1.2)

  loadings <- fit_draws(fit)$loadings
  expect_equal(dim(loadings), c(3, 2, posterior::ndraws(fit$draws("P"))))
  expect_length(draws_mean(fit, "factor_magnitude"), 2)
  expect_true(all(draws_mean(fit, "uniqueness") > 0))
})

test_that("factor P model validates its options", {
  sim <- simulate_animal(n_fam = 1, size = 100, E = diag(3), seed = 52,
                         genetic = FALSE)
  expect_error(
    lm_phen_cov(sim$Y, sim$X, covariance_structure = "factor",
                n_factors = 3),
    "`n_factors` must be smaller"
  )
  expect_error(
    lm_phen_cov(sim$Y, sim$X, covariance_structure = "factor",
                n_factors = 2, model_file = tempfile()),
    "`model_file` cannot be combined"
  )
  expect_error(
    lm_phen_cov(sim$Y, sim$X, n_factors = 2),
    "`n_factors` is only used"
  )
  expect_error(
    lm_phen_cov(sim$Y, sim$X, column_decay = 0.5),
    "`column_decay` is only used"
  )
})
