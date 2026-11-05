test_that("a custom model file fit cannot have its prior reconstructed", {
  skip_if_no_cmdstan()
  # the same program the built-in univariate model uses, so this costs no extra
  # compile: cstan() hashes the code and reuses the cached binary
  path <- system.file("stan", "animalModelUni.stan", package = "stanimal")
  sim <- simulate_animal(n_fam = 10, size = 4, G = matrix(0.5), E = matrix(1),
                         seed = 75)
  fit <- suppressWarnings(lmm_animal(sim$Y[, 1], sim$X, sim$A,
                                     model_file = path, iter = 300, warmup = 150,
                                     chains = 1, cores = 1, seed = 1,
                                     refresh = 0, show_messages = FALSE))
  expect_true(attr(fit, "custom_model"))
  expect_error(stanimal:::.prior_posterior(fit, NULL, NULL, 100, NULL),
               "custom Stan file")

  prior <- animal_prior_draws(1, 100)
  expect_no_error(stanimal:::.prior_posterior(fit, NULL, NULL, 100, prior))
})

test_that("a built-in fit is not flagged as custom", {
  skip_if_no_cmdstan()
  sim <- simulate_animal(n_fam = 10, size = 4, G = matrix(0.5), E = matrix(1),
                         seed = 76)
  fit <- suppressWarnings(lmm_animal(sim$Y[, 1], sim$X, sim$A,
                                     iter = 300, warmup = 150, chains = 1,
                                     cores = 1, seed = 1, refresh = 0,
                                     show_messages = FALSE))
  expect_false(attr(fit, "custom_model"))
  expect_no_error(stanimal:::.prior_posterior(fit, NULL, NULL, 100, NULL))
})

test_that("lm_phen_cov records its covariance configuration", {
  skip_if_no_cmdstan()
  sim <- simulate_animal(n_fam = 1, size = 120, E = diag(2), seed = 77,
                         genetic = FALSE)
  plain <- suppressWarnings(lm_phen_cov(sim$Y, sim$X, iter = 200, warmup = 100,
                                        chains = 1, cores = 1, seed = 1,
                                        refresh = 0, show_messages = FALSE))
  expect_equal(attr(plain, "covariance_structure"), "unstructured")
  expect_equal(attr(plain, "lkj_prior"), 2)
  expect_false(attr(plain, "custom_model"))
  expect_null(attr(plain, "n_factors"))

  factor_fit <- suppressWarnings(lm_phen_cov(
    sim$Y, sim$X, covariance_structure = "factor", n_factors = 1,
    column_decay = 0.9, iter = 200, warmup = 100, chains = 1, cores = 1,
    seed = 1, refresh = 0, show_messages = FALSE))
  expect_equal(attr(factor_fit, "covariance_structure"), "factor")
  expect_equal(attr(factor_fit, "n_factors"), 1)
  expect_equal(attr(factor_fit, "column_decay"), 0.9)
})
