skip_if_no_cmdstan()

shared_factor_truth <- function() {
  lambda <- matrix(c(
    0.80,  0.00,
    0.70,  0.40,
    0.50,  0.00,
    0.00,  0.70,
    0.00,  0.50,
    0.60, -0.30
  ), 6, 2, byrow = TRUE)
  h2 <- c(0.8, 0.4)
  psi_G <- c(0.30, 0.35, 0.30, 0.40, 0.35, 0.30)
  psi_E <- c(0.50, 0.45, 0.50, 0.45, 0.50, 0.45)
  G <- lambda %*% diag(h2) %*% t(lambda) + diag(psi_G)
  E <- lambda %*% diag(1 - h2) %*% t(lambda) + diag(psi_E)
  list(G = G, E = E, h2 = h2, lambda = lambda)
}

test_that("shared-factor model samples and exposes factor heritabilities", {
  truth <- shared_factor_truth()
  set.seed(41)
  n_fam <- 30L; size <- 8L; n <- n_fam * size; K <- nrow(truth$G)
  A <- make_A(n_fam, size)
  X <- cbind(Int = 1, x = stats::rnorm(n, 50, 10),
             z = stats::rbinom(n, 1, 0.5))
  B <- rbind(rep(10, K), rep(0.05, K), rep(1, K))
  Y <- X %*% B +
    t(chol(A)) %*% rmvn(n, truth$G) + rmvn(n, truth$E)
  sim <- list(Y = Y, X = X, A = A)
  fit <- suppressWarnings(lmm_animal(
    sim$Y, sim$X, sim$A,
    genetic_structure = "factor", residual_structure = "shared_factor",
    n_factors = 2, hs_global_scale = 2, hs_slab_scale = 2,
    column_decay = 0.9, iter = 1200, warmup = 600,
    chains = 2, cores = 2, seed = 1, refresh = 0, show_messages = FALSE
  ))
  expect_s3_class(fit, "CmdStanMCMC")
  expect_equal(attr(fit, "residual_structure"), "shared_factor")
  expect_equal(attr(fit, "n_factors"), 2)
  diagnostics <- fit$diagnostic_summary(quiet = TRUE)
  expect_equal(sum(diagnostics$num_divergent), 0)
  summaries <- fit$summary()
  covariance <- summaries[grepl("^(G|E|h2)\\[", summaries$variable), ]
  expect_lt(max(covariance$rhat, na.rm = TRUE), 1.05)

  G_hat <- post_matrix(fit, "G", 6)
  E_hat <- post_matrix(fit, "E", 6)
  expect_lt(max(abs(G_hat - truth$G)), 1.2)
  expect_lt(max(abs(E_hat - truth$E)), 1.2)

  h2_factor <- fit$summary("factor_heritability")$mean
  expect_length(h2_factor, 2)
  expect_true(all(h2_factor >= 0 & h2_factor <= 1))
  expect_true(all(is.finite(h2_factor)))
  expect_length(draws_mean(fit, "factor_heritability"), 2)
  expect_length(draws_mean(fit, "factor_magnitude"), 2)

  loadings <- fit_draws(fit)$loadings
  expect_equal(dim(loadings), c(6, 2, posterior::ndraws(fit$draws("G"))))
  uniqueness <- draws_mean(fit, "uniqueness")
  expect_length(uniqueness, 6)
  expect_true(all(uniqueness > 0))
  expect_equal(dim(fit$draws("log_lik"))[3], nrow(sim$Y))
})

test_that("shared_factor requires the factor genetic structure", {
  sim <- simulate_animal(n_fam = 10, size = 4, G = diag(2), E = diag(2),
                         seed = 42)
  expect_error(
    lmm_animal(sim$Y, sim$X, sim$A, residual_structure = "shared_factor"),
    "requires `genetic_structure"
  )
  expect_error(
    lmm_animal(sim$Y, sim$X, sim$A,
               genetic_structure = "unstructured",
               residual_structure = "shared_factor", n_factors = 2),
    "requires `genetic_structure"
  )
  expect_error(
    lmm_animal(sim$Y, sim$X, sim$A,
               genetic_structure = "factor",
               residual_structure = "shared_factor",
               n_factors = 1, model_file = tempfile()),
    "`model_file` cannot be combined"
  )
  expect_error(
    lmm_animal(sim$Y, sim$X, sim$A,
               genetic_structure = "factor",
               residual_structure = "shared_factor",
               n_factors = 1, beta_prior = c(1, -1)),
    "`beta_prior` must be two positive numbers"
  )
})
