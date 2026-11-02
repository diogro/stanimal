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

test_that("factor P describes the residuals, not the marginal covariance", {
  P <- phenotypic_factor_truth()
  sim <- simulate_animal(n_fam = 1, size = 400, E = P, seed = 53,
                         genetic = FALSE)
  # A fixed effect far larger than the residual variance, so the marginal
  # covariance of Y is much bigger than the residual covariance. If the model
  # pinned the diagonal of P to the marginal variance it could only fit the
  # likelihood by inflating the correlations.
  X <- sim$X
  B <- matrix(0, nrow(P), ncol(X))
  B[2L, ] <- c(6, -5, 4) * sqrt(diag(P))
  Y <- X %*% B + matrix(stats::rnorm(nrow(X) * 3), nrow(X)) %*% chol(P)
  colnames(Y) <- colnames(sim$Y)

  fit <- suppressWarnings(lm_phen_cov(
    Y, X, covariance_structure = "factor", n_factors = 2,
    hs_global_scale = 2, hs_slab_scale = 2, column_decay = 0.9,
    iter = 1000, warmup = 500, chains = 2, cores = 2, seed = 1,
    refresh = 0, show_messages = FALSE
  ))
  residual_sd <- sqrt(diag(P))
  marginal_sd <- apply(Y, 2L, stats::sd)
  sigma_P <- draws_mean(fit, "sigma_P")
  expect_lt(max(abs(sigma_P - residual_sd)), 0.25)
  expect_gt(min(marginal_sd / sigma_P), 2)
  expect_lt(max(abs(stats::cov2cor(draws_mean(fit, "P")) -
                      stats::cov2cor(P))), 0.2)
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
