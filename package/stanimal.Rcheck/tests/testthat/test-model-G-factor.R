skip_if_no_cmdstan()

fit_factor <- function(sim, iter = 1200, warmup = floor(iter / 2), ...) {
  suppressWarnings(lmm_animal(
    sim$Y, sim$X, sim$A, genetic_structure = "factor", n_factors = 2,
    hs_global_scale = 2, hs_slab_scale = 2, iter = iter, warmup = warmup,
    chains = 2, cores = 2, seed = 1, refresh = 0, show_messages = FALSE, ...
  ))
}

true_factor_G <- function() {
  lambda <- matrix(c(0.8, 0, 0.7, 0.3, 0.6, -0.4), 3, 2, byrow = TRUE)
  lambda %*% t(lambda) + diag(c(0.36, 0.42, 0.48))
}

test_that("factor G model samples and exposes factor summaries", {
  G <- true_factor_G()
  sim <- simulate_animal(n_fam = 30, size = 8, G = G,
                         E = diag(c(0.8, 0.9, 0.7)), seed = 31)
  colnames(sim$Y) <- c("t1", "t2", "t3")
  fit <- fit_factor(sim)
  expect_s3_class(fit, "CmdStanMCMC")
  expect_equal(attr(fit, "genetic_structure"), "factor")
  expect_equal(attr(fit, "n_factors"), 2)

  expect_clean_fit(fit)
  # the posterior mean recovers G to about 0.4 on entries of order one
  expect_lt(max(abs(post_matrix(fit, "G", 3) - G)), 0.6)

  loadings <- fit_draws(fit)$loadings
  expect_equal(dim(loadings), c(3, 2, posterior::ndraws(fit$draws("G"))))
  expect_equal(dimnames(loadings)[[1]], colnames(sim$Y))
  expect_equal(dimnames(loadings)[[2]], c("factor1", "factor2"))
  expect_equal(dim(draws_mean(fit, "loadings")), c(3, 2))

  uniqueness <- draws_mean(fit, "uniqueness")
  communality <- draws_mean(fit, "communality")
  expect_true(all(uniqueness > 0 & uniqueness < 1))
  expect_true(all(communality > 0 & communality < 1))
  expect_equal(unname(uniqueness + communality), rep(1, 3), tolerance = 1e-6)

  s_G <- draws_mean(fit, "s_G")
  expect_true(all(s_G > 0))
  expect_true(all(is.finite(fit$summary("loadings")$mean)))
  ll <- fit$draws("log_lik")
  expect_equal(dim(ll)[3], nrow(sim$Y))
  expect_true(all(is.finite(ll)))

  # The model is parameterized by h2 and a phenotypic scale, so the reported h2
  # is a parameter, and it is still exactly the diagonal ratio.
  dm <- as.matrix(fit$draws(c("G", "P", "h2", "s_P"),
                            format = "draws_matrix"))
  expect_equal(unname(dm[, sprintf("h2[%d]", 1:3)]),
               unname(dm[, sprintf("G[%d,%d]", 1:3, 1:3)] /
                        dm[, sprintf("P[%d,%d]", 1:3, 1:3)]),
               tolerance = 1e-6)
  expect_true(all(dm[, sprintf("s_P[%d]", 1:3)] > 0))
  expect_true(all(dm[, sprintf("h2[%d]", 1:3)] > 0 &
                    dm[, sprintf("h2[%d]", 1:3)] < 1))
  expect_equal(attr(fit, "parametrization"), "h2")
  expect_equal(attr(fit, "beta_prior"), c(1, 1))
})

test_that("factor model uses the sublinear n_factors default", {
  expect_equal(stanimal:::.default_n_factors(2), 1)
  expect_equal(stanimal:::.default_n_factors(3), 2)
  expect_equal(stanimal:::.default_n_factors(5), 2)
  expect_equal(stanimal:::.default_n_factors(10), 3)
  expect_equal(stanimal:::.default_n_factors(100), 10)
  expect_equal(stanimal:::.default_n_factors(1000), 20)

  sim <- simulate_animal(n_fam = 10, size = 4, G = true_factor_G(),
                         E = diag(3), seed = 32)
  expect_warning(
    fit_default <- lmm_animal(
      sim$Y, sim$X, sim$A, genetic_structure = "factor",
      iter = 200, warmup = 100, chains = 1, cores = 1,
      seed = 1, refresh = 0, show_messages = FALSE),
    "n_factors"
  )
  expect_equal(attr(fit_default, "n_factors"), 2)
})

test_that("factor model validates n_factors and model files", {
  sim <- simulate_animal(n_fam = 10, size = 4, G = true_factor_G(),
                         E = diag(3), seed = 32)
  expect_error(
    lmm_animal(sim$Y, sim$X, sim$A, genetic_structure = "factor",
               n_factors = 3),
    "`n_factors` must be smaller"
  )
  expect_error(
    lmm_animal(sim$Y, sim$X, sim$A, genetic_structure = "factor",
               n_factors = 2, model_file = tempfile()),
    "`model_file` cannot be combined"
  )
  expect_warning(
    lmm_animal(sim$Y, sim$X, sim$A, column_decay = 0.5),
    "`column_decay` is ignored"
  )
  expect_error(
    lmm_animal(sim$Y, sim$X, sim$A, genetic_structure = "factor",
               n_factors = 2, column_decay = 0),
    "column_decay"
  )
  expect_error(
    lmm_animal(sim$Y, sim$X, sim$A, genetic_structure = "factor",
               n_factors = 2, column_decay = 1.2),
    "column_decay"
  )
  expect_error(
    lmm_animal(sim$Y[, 1], sim$X, sim$A, genetic_structure = "factor",
               n_factors = 2),
    "at least two traits"
  )
})

test_that("column_decay is passed to the factor model", {
  sim <- simulate_animal(n_fam = 20, size = 6, G = true_factor_G(),
                         E = diag(c(0.8, 0.9, 0.7)), seed = 34)
  fit <- fit_factor(sim, iter = 400, column_decay = 0.5)
  expect_s3_class(fit, "CmdStanMCMC")
  magnitude <- fit$summary("factor_magnitude")$mean
  expect_length(magnitude, 2)
  expect_true(all(is.finite(magnitude)))
  expect_true(all(is.finite(fit$summary("loadings")$mean)))
})

test_that("factor model handles a diagonal residual structure", {
  sim <- simulate_animal(n_fam = 20, size = 6, G = true_factor_G(),
                         E = diag(c(0.8, 0.9, 0.7)), seed = 33)
  fit <- fit_factor(sim, iter = 600, residual_structure = "diagonal")
  expect_s3_class(fit, "CmdStanMCMC")
  E <- post_matrix(fit, "E", 3)
  expect_equal(E[upper.tri(E)], rep(0, 3), tolerance = 1e-8)
  expect_true(all(fit$summary("uniqueness")$mean > 0))
})

test_that("the factor model puts beta_prior on the trait heritabilities", {
  sim <- simulate_animal(n_fam = 10, size = 5, G = true_factor_G(),
                         E = diag(c(0.8, 0.9, 0.7)), seed = 35)
  short <- function(prior) {
    lmm_animal(sim$Y, sim$X, sim$A, genetic_structure = "factor",
               n_factors = 2, hs_global_scale = 2, beta_prior = prior,
               iter = 200, warmup = 100, chains = 1, cores = 1, seed = 1,
               refresh = 0, show_messages = FALSE)
  }
  # the prior is used by this model, so it is not reported as ignored
  warns <- character()
  fit_flat <- withCallingHandlers(
    short(c(1, 1)),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_false(any(grepl("ignored", warns)))
  expect_equal(attr(fit_flat, "beta_prior"), c(1, 1))

  # Beta(5, 50) has mean 0.091, so a live prior must drag h2 down
  fit_low <- suppressWarnings(short(c(5, 50)))
  expect_true(all(draws_mean(fit_low, "h2") < draws_mean(fit_flat, "h2")))
  expect_true(all(draws_mean(fit_low, "h2") < 0.4))
})
