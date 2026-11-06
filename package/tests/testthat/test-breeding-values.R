# Estimated breeding values.
#
# The animal models integrate the additive genetic effects out, so `u` is never
# a sampled quantity. `breeding_values()` recovers it by solving the conditional
# mean at the posterior means of G, E and beta. The first three checks verify
# that solve against the definition, written out with Kronecker products, which
# is the only way to be sure the eigendecomposition shortcut is the same
# quantity.

test_that("the predicted values are the conditional means of u", {
  set.seed(11)
  n <- 8
  k <- 3
  A <- make_A(2, 4)
  X <- cbind(Int = 1, x = stats::rnorm(n))
  beta <- rbind(c(2, 0.5), c(1, -0.2), c(0, 0.3))
  G <- matrix(c(1.0, 0.4, 0.1,
                0.4, 0.8, 0.2,
                0.1, 0.2, 0.6), 3)
  E <- diag(c(0.5, 0.4, 0.3))
  # kronecker(A, G) orders the traits of one individual together, which is a
  # traits by individuals vector; the response is the other way round, so the
  # draws are reshaped and transposed, and the reference below is kronecker(G, A)
  u <- t(matrix(rmvn(1, kronecker(A, G)), k, n))
  e <- t(matrix(rmvn(1, kronecker(E, diag(n))), k, n))
  Y <- X %*% t(beta) + u + e

  out <- stanimal:::.conditional_breeding_values(Y, X, A, G, E, beta)

  # as.vector() of an individuals by traits matrix puts trait 1 for every
  # individual first, so the covariance of that vector is kronecker(G, A)
  V <- kronecker(G, A) + kronecker(E, diag(n))
  Ag <- kronecker(G, A)
  residual <- as.vector(Y - X %*% t(beta))
  expected <- matrix(Ag %*% solve(V, residual), n, k)
  expect_equal(unname(out$values), expected, tolerance = 1e-10)

  # Var(u | y) = Cov(u, y) Var(y)^-1 Cov(y, u)
  conditional <- Ag - Ag %*% solve(V, Ag)
  expect_equal(unname(out$se), sqrt(matrix(diag(conditional), n, k)),
               tolerance = 1e-10)
  reliability <- matrix(1 - diag(conditional) /
                          (rep(diag(G), each = n) * diag(A)), n, k)
  expect_equal(unname(out$reliability), reliability, tolerance = 1e-10)
})

test_that("reliability is relative to the prior variance A[n, n] * G[k, k]", {
  G <- matrix(c(1, 0.3, 0.3, 1), 2)
  E <- diag(2)
  Y <- matrix(c(0.5, -0.2, 0.1, 0.4, -0.3, 0.2), 3)
  X <- matrix(1, 3, 1)
  beta <- c(0, 0)
  # with A = a I each individual is independent: Var(u | y) = a G - a G (a G + E)^-1 a G
  for (a in c(0.4, 2)) {
    out <- stanimal:::.conditional_breeding_values(Y, X, a * diag(3), G, E, beta)
    conditional <- a * G - a * G %*% solve(a * G + E, a * G)
    expected <- 1 - matrix(diag(conditional), 3, 2, byrow = TRUE) / (a * matrix(diag(G), 3, 2, byrow = TRUE))
    expect_equal(unname(out$reliability), expected, tolerance = 1e-10)
  }
})

test_that("with A = I the prediction is ridge shrinkage of each residual", {
  set.seed(12)
  n <- 6
  k <- 2
  A <- diag(n)
  X <- matrix(1, n, 1L)
  beta <- matrix(c(3, -1), k, 1L)
  G <- matrix(c(1.2, 0.3, 0.3, 0.9), 2)
  E <- matrix(c(0.8, 0.1, 0.1, 0.7), 2)
  u <- t(matrix(rmvn(1, kronecker(A, G)), k, n))
  e <- t(matrix(rmvn(1, kronecker(E, diag(n))), k, n))
  Y <- X %*% t(beta) + u + e

  out <- stanimal:::.conditional_breeding_values(Y, X, A, G, E, beta)

  # no relatives, so the only information about an individual is its own
  # residual, shrunk by G (G + E)^-1
  shrink <- G %*% solve(G + E)
  expected <- t(shrink %*% t(Y - X %*% t(beta)))
  expect_equal(unname(out$values), expected, tolerance = 1e-10)
  # and a conditional variance below the unconditional one
  expect_true(all(out$reliability > 0 & out$reliability < 1))
  expect_true(all(out$se <= sqrt(rep(diag(G), each = n)) + 1e-10))
})

test_that("a prediction keeps the individual and trait labels", {
  set.seed(13)
  n <- 5
  k <- 2
  A <- diag(n)
  X <- matrix(1, n, 1L)
  Y <- matrix(stats::rnorm(n * k), n, k,
              dimnames = list(paste0("id", seq_len(n)), c("a", "b")))
  out <- stanimal:::.conditional_breeding_values(Y, X, A, diag(k), diag(k),
                                                matrix(0, k, 1L))
  expect_equal(rownames(out$values), paste0("id", seq_len(n)))
  expect_equal(colnames(out$values), c("a", "b"))
  expect_equal(dim(out$se), c(n, k))
  expect_equal(dim(out$reliability), c(n, k))
})

test_that("only an animal fit has breeding values", {
  n <- 6
  Y <- matrix(stats::rnorm(n * 3), n, 3)
  X <- matrix(1, n, 1L)
  phenotypic <- structure(list(),
                          stanimal_data = list(family = "phenotypic", Y = Y,
                                               X = X, A = NULL))
  expect_error(breeding_values(phenotypic), "needs an animal model")

  expect_error(breeding_values(list(a = 1)), "does not carry the data")

  # the Gibbs sampler does not record the data it was fitted to
  gibbs <- lmm_animal_gibbs(Y, X, make_A(2, 3), iter = 40, warmup = 20)
  expect_error(breeding_values(gibbs), "does not carry the data")
})

test_that("univariate plug-in variances are the means of squared SD draws", {
  m <- cbind("sigma_G" = c(1, 3), "sigma_E" = c(2, 4), "beta[1]" = c(0, 0))
  Y <- matrix(c(-1, 1), 2, 1, dimnames = list(c("a", "b"), "trait"))
  fit <- structure(list(draws = function(...) posterior::as_draws_matrix(m)),
                   class = "CmdStanMCMC", y_names = "trait", x_names = "Int",
                   stanimal_data = list(family = "animal", Y = Y,
                                        X = matrix(1, 2, 1), A = diag(2)))
  ebv <- breeding_values(fit)
  expect_equal(unname(ebv$G), matrix(5))
  expect_equal(unname(ebv$E), matrix(10))
  expect_equal(unname(ebv$values), Y / 3, ignore_attr = TRUE)
  expect_equal(unname(ebv$se), matrix(sqrt(10 / 3), 2, 1))
  expect_equal(unname(ebv$reliability), matrix(1 / 3, 2, 1))
})

test_that("missing-response uncertainty conditions on observed entries", {
  A <- matrix(c(1.3, 0.4, 0.2, 0.4, 1, 0.3, 0.2, 0.3, 0.8), 3)
  G <- matrix(c(1.2, 0.4, 0.4, 0.9), 2)
  E <- matrix(c(0.7, -0.15, -0.15, 0.6), 2)
  Y <- matrix(c(0.3, -0.4, 0.6, 1.1, -0.2, 0.7), 3, 2)
  X <- matrix(1, 3, 1)
  beta <- matrix(c(0.2, -0.1), 2, 1)
  K <- kronecker(G, A)
  V <- K + kronecker(E, diag(3))
  for (missing in list(matrix(c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), 3),
                       matrix(c(TRUE, FALSE, FALSE, FALSE, TRUE, FALSE), 3),
                       matrix(c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE), 3),
                       matrix(c(TRUE, FALSE, TRUE, FALSE, TRUE, TRUE), 3))) {
    observed <- which(!as.vector(missing))
    expected <- K - K[, observed, drop = FALSE] %*%
      solve(V[observed, observed, drop = FALSE], K[observed, , drop = FALSE])
    out <- stanimal:::.conditional_breeding_values(
      Y, X, A, G, E, beta, missing = missing)
    complete <- stanimal:::.conditional_breeding_values(Y, X, A, G, E, beta)
    expect_equal(out$values, complete$values)
    expect_equal(unname(out$se), matrix(sqrt(diag(expected)), 3, 2),
                 tolerance = 1e-10)
    expect_equal(unname(out$reliability),
                 matrix(1 - diag(expected) / diag(K), 3, 2), tolerance = 1e-10)
    expect_true(all(out$se >= complete$se - 1e-10))

    mu <- as.vector(X %*% t(beta))
    filled <- mu + V[, observed, drop = FALSE] %*%
      solve(V[observed, observed, drop = FALSE],
            as.vector(Y)[observed] - mu[observed])
    expected_mean <- K[, observed, drop = FALSE] %*%
      solve(V[observed, observed, drop = FALSE],
            as.vector(Y)[observed] - mu[observed])
    imputed <- stanimal:::.conditional_breeding_values(
      matrix(filled, 3, 2), X, A, G, E, beta, missing = missing)
    expect_equal(unname(imputed$values), matrix(expected_mean, 3, 2),
                 tolerance = 1e-10)
  }
  all_missing <- matrix(TRUE, 3, 2)
  out <- stanimal:::.conditional_breeding_values(
    Y, X, A, G, E, beta, missing = all_missing)
  expect_equal(unname(out$se), matrix(sqrt(diag(K)), 3, 2), tolerance = 1e-10)
  expect_equal(unname(out$reliability), matrix(0, 3, 2), tolerance = 1e-10)
})

test_that("an unrelated individual with no observations has zero reliability", {
  Y <- matrix(c(0, 1, -1), 3, 1)
  missing <- matrix(c(TRUE, FALSE, FALSE), 3, 1)
  out <- stanimal:::.conditional_breeding_values(
    Y, matrix(1, 3, 1), diag(3), matrix(1), matrix(1), matrix(0),
    missing = missing)
  expect_equal(out$se[1, 1], 1)
  expect_equal(out$reliability[1, 1], 0, tolerance = 1e-12)
  expect_equal(out$se[2, 1], sqrt(0.5))
  expect_equal(out$reliability[2, 1], 0.5)
})

test_that("breeding_values uses imputations for means and the observed mask for SE", {
  Y <- matrix(c(NA, -1, 1), 3, 1)
  m <- cbind("sigma_G" = c(1, 1), "sigma_E" = c(1, 1),
             "beta[1]" = c(0, 0))
  fit <- structure(list(draws = function(...) posterior::as_draws_matrix(m)),
                   class = "CmdStanMCMC", y_names = "trait", x_names = "Int",
                   stanimal_data = list(family = "animal", Y = Y,
                                        X = matrix(1, 3, 1), A = diag(3)))
  local_mocked_bindings(
    imputed_data = function(fit, ndraws, seed) {
      expect_equal(ndraws, 2)
      expect_equal(seed, 17)
      list(Y = array(c(2, -1, 1, 4, -1, 1), c(3, 1, 2)))
    }, .package = "stanimal")
  out <- breeding_values(fit, ndraws = 2, seed = 17)
  expect_equal(out$values[1, 1], 1.5)
  expect_equal(out$se[1, 1], 1)
  expect_equal(out$reliability[1, 1], 0, tolerance = 1e-12)
  no_se <- breeding_values(fit, se = FALSE, ndraws = 2, seed = 17)
  expect_equal(no_se$values, out$values)
  expect_null(no_se$se)
  expect_null(no_se$reliability)
})

test_that("se FALSE does not compute the missing-data covariance correction", {
  local_mocked_bindings(
    .missing_breeding_variance = function(...) stop("unexpected correction"),
    .package = "stanimal")
  expect_silent(stanimal:::.conditional_breeding_values(
    matrix(c(1, 2), 2), matrix(1, 2, 1), diag(2), matrix(1), matrix(1),
    matrix(0), se = FALSE, missing = matrix(c(TRUE, FALSE), 2)))
})

skip_if_no_cmdstan()

ebv_fit <- function(sim, iter = 800, ...) {
  suppressWarnings(lmm_animal(sim$Y, sim$X, sim$A, iter = iter,
                              warmup = floor(iter / 2), chains = 2, cores = 2,
                              seed = 5, refresh = 0, show_messages = FALSE, ...))
}

test_that("breeding values from a fit track the simulated ones", {
  sim <- simulate_animal(n_fam = 25, size = 8, G = G3, E = E3, seed = 31)
  fit <- ebv_fit(sim)
  ebv <- breeding_values(fit)

  expect_s3_class(ebv, "stanimal_breeding_values")
  expect_equal(dim(ebv$values), dim(sim$Y))
  expect_equal(colnames(ebv$values), paste0("trait", 1:3))
  expect_true(all(ebv$reliability >= 0 & ebv$reliability <= 1))
  expect_true(all(ebv$se > 0))
  # no prediction can be more uncertain than the trait's additive variance
  expect_true(all(ebv$se <= sqrt(rep(diag(G3), each = nrow(sim$Y))) + 1e-8))

  # the prediction is correlated with the simulated breeding value, and closer
  # to it than the individual's own deviation from the fixed effects
  residual <- sim$Y - sim$X %*% t(draws_mean(fit, "beta"))
  rmse <- function(x) sqrt(mean(x^2))
  for (k in seq_len(ncol(sim$Y))) {
    expect_gt(stats::cor(ebv$values[, k], sim$u[, k]), 0.6)
    expect_lt(rmse(ebv$values[, k] - sim$u[, k]),
              rmse(residual[, k] - sim$u[, k]))
  }
})

test_that("a univariate fit gives breeding values too", {
  sim <- simulate_animal(n_fam = 20, size = 6, G = G3, E = E3, seed = 32)
  fit <- ebv_fit(sim, iter = 600)
  fit1 <- suppressWarnings(lmm_animal(sim$Y[, 1], sim$X, sim$A, iter = 600,
                                      warmup = 300, chains = 2, cores = 2,
                                      seed = 5, refresh = 0,
                                      show_messages = FALSE))
  ebv <- breeding_values(fit1)
  expect_equal(dim(ebv$values), c(nrow(sim$Y), 1L))
  expect_true(all(ebv$reliability >= 0 & ebv$reliability <= 1))
  # the univariate model has less information, so the predictions shrink more
  multi <- breeding_values(fit)
  expect_true(stats::sd(ebv$values[, 1]) <= stats::sd(multi$values[, 1]) + 0.2)
})

test_that("factor-structured fits predict from the G they report", {
  sim <- simulate_animal(n_fam = 15, size = 6, G = G3, E = E3, seed = 35)
  factor_fit <- ebv_fit(sim, iter = 600, genetic_structure = "factor",
                        n_factors = 1)
  ebv <- breeding_values(factor_fit)
  expect_equal(dim(ebv$values), dim(sim$Y))
  expect_true(all(ebv$reliability >= 0 & ebv$reliability <= 1))
  # a factor fit reports G on the response scale like any other, so the pass
  # uses it directly instead of needing the loadings
  expect_true(all(ebv$se <= sqrt(rep(diag(G3), each = nrow(sim$Y))) + 1e-8))

  shared_fit <- ebv_fit(sim, iter = 600, genetic_structure = "factor",
                        residual_structure = "shared_factor", n_factors = 1)
  shared <- breeding_values(shared_fit)
  expect_equal(dim(shared$values), dim(sim$Y))
  expect_true(all(shared$se > 0))
})

test_that("missing values are predicted through their conditional mean", {
  sim <- simulate_animal(n_fam = 20, size = 6, G = G3, E = E3, seed = 33)
  incomplete <- sim$Y
  incomplete[c(3, 17, 25, 41), 1] <- NA
  incomplete[c(8, 30), 3] <- NA
  fit <- ebv_fit(list(Y = incomplete, X = sim$X, A = sim$A), iter = 600)
  ebv <- breeding_values(fit, ndraws = 200, seed = 4)

  expect_equal(dim(ebv$values), dim(sim$Y))
  expect_false(anyNA(ebv$values))
  expect_true(all(ebv$se > 0))
  # a mouse whose trait was not measured is still predicted, from its relatives
  expect_false(any(ebv$values[c(3, 17), 1] == 0))
  # and its reliability is below that of a mouse with its own record
  expect_lt(mean(ebv$reliability[c(3, 17, 25, 41), 1]),
            mean(ebv$reliability[-c(3, 17, 25, 41), 1]))
  completed <- stanimal:::.completed_mean(imputed_data(fit, ndraws = 200, seed = 4))
  baseline <- stanimal:::.conditional_breeding_values(
    completed, sim$X, sim$A, ebv$G, ebv$E, draws_mean(fit, "beta"))
  expect_equal(ebv$values, baseline$values)
  K <- kronecker(ebv$G, sim$A)
  V <- K + kronecker(ebv$E, diag(nrow(incomplete)))
  observed <- which(!is.na(as.vector(incomplete)))
  covariance <- K - K[, observed, drop = FALSE] %*%
    solve(V[observed, observed, drop = FALSE], K[observed, , drop = FALSE])
  expect_equal(unname(ebv$se),
               matrix(sqrt(diag(covariance)), nrow(incomplete)),
               tolerance = 1e-8)
  expect_equal(unname(ebv$reliability),
               matrix(1 - diag(covariance) / diag(K), nrow(incomplete)),
               tolerance = 1e-8)
})

test_that("se = FALSE skips the uncertainty", {
  sim <- simulate_animal(n_fam = 12, size = 5, G = G3, E = E3, seed = 34)
  fit <- ebv_fit(sim, iter = 500)
  ebv <- breeding_values(fit, se = FALSE)
  expect_null(ebv$se)
  expect_null(ebv$reliability)
  expect_equal(ebv$values, breeding_values(fit)$values)
})
