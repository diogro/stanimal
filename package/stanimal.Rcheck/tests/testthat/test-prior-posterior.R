# The prior drawn by plot_prior_posterior() has to be the prior the fit used,
# not always the unstructured one.

prior_only <- function(fit, ...) {
  d <- stanimal:::.prior_posterior(fit, NULL, NULL, 500, NULL)
  d$value[d$source == "Prior"]
}

test_that("a fit records the configuration its prior needs", {
  skip_if_no_cmdstan()
  sim <- simulate_animal(n_fam = 10, size = 4, E = diag(2), seed = 71,
                         genetic = FALSE)
  fit <- suppressWarnings(lm_phen_cov(sim$Y, sim$X, iter = 300, warmup = 150,
                                      chains = 1, cores = 1, seed = 1,
                                      refresh = 0, show_messages = FALSE))
  expect_equal(attr(fit, "y_names"), colnames(sim$Y))
})

test_that("the h2 parametrization is compared with its beta prior", {
  skip_if_no_cmdstan()
  sim <- simulate_animal(n_fam = 10, size = 4, G = matrix(0.5), E = matrix(1),
                         seed = 72)
  fit <- suppressWarnings(lmm_animal(sim$Y[, 1], sim$X, sim$A,
                                     parametrization = "h2",
                                     beta_prior = c(2, 8),
                                     iter = 300, warmup = 150, chains = 1,
                                     cores = 1, seed = 1, refresh = 0,
                                     show_messages = FALSE))
  expect_equal(attr(fit, "parametrization"), "h2")
  expect_equal(attr(fit, "beta_prior"), c(2, 8))

  set.seed(1)
  draws <- prior_only(fit)
  # beta(2, 8) has mean 0.2; the scaled half-normal prior would be near 0.5
  expect_lt(abs(mean(draws) - 0.2), 0.03)
  expect_true(all(draws > 0 & draws < 1))
})

test_that("an unstructured fit is compared with the LKJ shape it was given", {
  skip_if_no_cmdstan()
  sim <- simulate_animal(n_fam = 10, size = 4, G = diag(2), E = diag(2),
                         seed = 73)
  fits <- lapply(c(1, 6), function(eta)
    suppressWarnings(lmm_animal(sim$Y, sim$X, sim$A, lkj_prior = eta,
                                iter = 300, warmup = 150, chains = 1, cores = 1,
                                seed = 1, refresh = 0, show_messages = FALSE)))
  expect_equal(vapply(fits, function(f) attr(f, "lkj_prior"), 0), c(1, 6))

  # the prior of the second fit is LKJ(6): its correlations are much tighter
  set.seed(1); loose <- prior_only(fits[[1]])
  set.seed(1); tight <- prior_only(fits[[2]])
  expect_lt(sd(tight), sd(loose))
})

test_that("a factor fit has no closed-form prior to compare against", {
  skip_if_no_cmdstan()
  sim <- simulate_animal(n_fam = 10, size = 4, G = diag(c(1, 0.5)),
                         E = diag(2), seed = 74)
  fit <- suppressWarnings(lmm_animal(sim$Y, sim$X, sim$A,
                                     genetic_structure = "factor",
                                     n_factors = 1, iter = 300, warmup = 150,
                                     chains = 1, cores = 1, seed = 1,
                                     refresh = 0, show_messages = FALSE))
  expect_error(stanimal:::.prior_posterior(fit, NULL, NULL, 100, NULL),
               "sparse factor model")
  # an explicit prior is still accepted, so the plot stays available
  set.seed(1)
  prior <- animal_prior_draws(2, 100, trait_names = colnames(sim$Y))
  expect_no_error(stanimal:::.prior_posterior(fit, NULL, NULL, 100, prior))
})

test_that("the beta prior applies to every trait", {
  set.seed(1)
  d2 <- animal_prior_draws(2, n = 50, beta_prior = c(1, 1))
  expect_setequal(unique(d2$type), c("Heritability", "Genetic correlation"))
  expect_true(all(d2$value[d2$type == "Heritability"] > 0 &
                    d2$value[d2$type == "Heritability"] < 1))
  # the correlations are still the LKJ, so they are in [-1, 1] and not zero
  rho <- d2$value[d2$type == "Genetic correlation"]
  expect_true(all(abs(rho) <= 1))
  expect_gt(stats::sd(rho), 0)
  expect_error(animal_prior_draws(1, n = 10, beta_prior = c(1, -1)),
               "`beta_prior`")
  d <- animal_prior_draws(1, n = 50, beta_prior = c(1, 1))
  expect_equal(unique(d$type), "Heritability")
  expect_true(all(d$value > 0 & d$value < 1))
  # a tight beta is reproduced, which is what makes the plot a prior check
  dt <- animal_prior_draws(1, n = 4000, beta_prior = c(20, 80))
  expect_equal(mean(dt$value), 0.2, tolerance = 0.02)
})
