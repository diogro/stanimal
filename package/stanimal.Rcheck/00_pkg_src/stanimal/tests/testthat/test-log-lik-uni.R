skip_if_no_cmdstan()

# The univariate models are the only ones without a pointwise log likelihood, so
# LOO and WAIC could not be computed for the fits people run most.
uni_log_lik <- function(fit) {
  draws <- fit$draws("log_lik", format = "draws_matrix")
  as.matrix(draws)
}

test_that("univariate fits report a pointwise log likelihood", {
  sim <- simulate_animal(n_fam = 20, size = 8, G = matrix(0.6), E = matrix(0.6),
                         seed = 91)
  fit <- suppressWarnings(lmm_animal(sim$Y, sim$X, sim$A, iter = 600,
                                     chains = 2, cores = 2, seed = 1,
                                     refresh = 0, show_messages = FALSE))
  ll <- uni_log_lik(fit)
  expect_equal(ncol(ll), nrow(sim$Y))
  expect_true(all(is.finite(ll)))
  # a per-individual log density of a plausible magnitude, not a placeholder
  expect_true(all(colMeans(ll) < 0))
})

test_that("the h2 parametrization reports the same density", {
  # both parametrizations describe the same model, so the total log likelihood
  # evaluated at the posterior mean has to agree; a tolerance of a few percent
  # leaves room for the two samplers visiting slightly different regions
  sim <- simulate_animal(n_fam = 20, size = 8, G = matrix(0.6), E = matrix(0.6),
                         seed = 92)
  variance_fit <- suppressWarnings(lmm_animal(
    sim$Y, sim$X, sim$A, parametrization = "variance", iter = 800, chains = 2,
    cores = 2, seed = 1, refresh = 0, show_messages = FALSE))
  h2_fit <- suppressWarnings(lmm_animal(
    sim$Y, sim$X, sim$A, parametrization = "h2", beta_prior = c(2, 2),
    iter = 800, chains = 2, cores = 2, seed = 1, refresh = 0,
    show_messages = FALSE))

  total_variance <- mean(rowSums(uni_log_lik(variance_fit)))
  total_h2 <- mean(rowSums(uni_log_lik(h2_fit)))
  expect_lt(abs(total_variance - total_h2) / abs(total_variance), 0.05)
})

test_that("the log likelihood matches the density the model block uses", {
  # A = I removes the pedigree, so lambda = 1 for everyone, the rotation is the
  # identity, and the marginal variance is sigma_G^2 + sigma_E^2 for every
  # individual. That makes the expected log likelihood computable from the fit's
  # own posterior means of the standard deviations alone.
  set.seed(93)
  n <- 200
  x <- stats::rnorm(n)
  X <- cbind(Int = 1, x = x)
  y <- 10 + 2 * x + stats::rnorm(n, sd = sqrt(0.5)) + stats::rnorm(n, sd = sqrt(0.5))
  fit <- suppressWarnings(lmm_animal(y, X, diag(n), iter = 800, chains = 2,
                                     cores = 2, seed = 1, refresh = 0,
                                     show_messages = FALSE))
  ll <- uni_log_lik(fit)

  sigma_G <- mean(fit$draws("sigma_G", format = "draws_matrix"))
  sigma_E <- mean(fit$draws("sigma_E", format = "draws_matrix"))
  beta <- fit$summary("beta", "mean")$mean
  y_sd <- stats::sd(y)
  residuals_standardized <- (y - X %*% beta) / y_sd
  expected <- stats::dnorm(residuals_standardized, 0,
                           sqrt(sigma_G^2 + sigma_E^2) / y_sd, log = TRUE)
  # compare the totals: a per-individual match would need the standardized-scale
  # parameters, which the model does not expose
  total_expected <- sum(expected)
  expect_lt(abs(mean(rowSums(ll)) - total_expected) / abs(total_expected), 0.05)
})
