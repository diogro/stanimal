# The prior drawn by plot_prior_posterior() has to be the prior the fit used,
# not always the unstructured one.

prior_only <- function(fit, ...) {
  d <- stanimal:::.prior_posterior(fit, NULL, NULL, 500, NULL)
  d$value[d$source == "Prior"]
}

prior_test_fit <- function(k = 2L, parametrization = "variance",
                           residual = "unstructured") {
  vars <- c(sprintf("h2[%d]", seq_len(k)),
            if (k > 1L) "corrG[1,2]")
  values <- matrix(rep(seq_along(vars) / 10, each = 20), 20,
                   dimnames = list(NULL, vars))
  fit <- posterior::as_draws_matrix(values)
  attr(fit, "genetic_structure") <- "unstructured"
  attr(fit, "residual_structure") <- residual
  attr(fit, "parametrization") <- parametrization
  attr(fit, "custom_model") <- FALSE
  attr(fit, "lkj_prior") <- 3
  attr(fit, "beta_prior") <- c(2, 8)
  attr(fit, "y_names") <- paste0("trait", seq_len(k))
  fit
}

test_that("automatic priors follow each unstructured model configuration", {
  for (k in 1:2) for (param in c("variance", "h2")) {
    for (residual in c("unstructured", "diagonal")) {
      fit <- prior_test_fit(k, param, residual)
      set.seed(104)
      d <- stanimal:::.prior_posterior(fit, NULL, NULL, 4000, NULL)
      prior <- d[d$source == "Prior", ]
      h2 <- prior$value[prior$type == "Heritability"]
      expect_lt(abs(mean(h2) - if (param == "h2") 0.2 else 0.5), 0.025)
      expect_lt(abs(var(h2) - if (param == "h2") 0.2 * 0.8 / 11 else 0.125),
                0.015)
      if (k == 2L) {
        rho <- prior$value[prior$type == "Genetic correlation"]
        expect_lt(abs(var(rho) - 1 / 7), 0.02)
      }
    }
  }
})

test_that("unknown or incomplete prior metadata never selects a default", {
  fit <- prior_test_fit()
  for (field in c("genetic_structure", "residual_structure", "parametrization",
                  "custom_model", "lkj_prior")) {
    incomplete <- fit
    attr(incomplete, field) <- NULL
    expect_error(plot_prior_posterior(incomplete, n_prior = 20),
                 "cannot.*reconstruct|cannot be reconstructed")
    expect_error(plot_forest(incomplete, n_prior = 20),
                 "cannot.*reconstruct|cannot be reconstructed")
  }
  h2fit <- prior_test_fit(parametrization = "h2")
  attr(h2fit, "beta_prior") <- NULL
  expect_error(plot_prior_posterior(h2fit, n_prior = 20), "cannot.*reconstruct")
  bare <- posterior::as_draws_matrix(matrix(as.numeric(fit), nrow(fit),
                                           dimnames = list(NULL, colnames(fit))))
  expect_error(plot_prior_posterior(bare, n_prior = 20), "cannot.*reconstruct")
  expect_error(plot_prior_posterior(bare, n_prior = 20, lkj_prior = 3),
               "cannot.*reconstruct")
  prior <- animal_prior_draws(2, 50)
  expect_s3_class(plot_prior_posterior(bare, prior = prior), "ggplot")
})

test_that("unsupported fits require explicit prior draws", {
  fit <- prior_test_fit()
  for (field in c("genetic_structure", "residual_structure", "custom_model")) {
    other <- fit
    attr(other, field) <- switch(field, genetic_structure = "factor",
                                 residual_structure = "shared_factor",
                                 custom_model = TRUE)
    expect_error(plot_prior_posterior(other, n_prior = 20), "prior")
    prior <- animal_prior_draws(2, 50, trait_names = c("trait1", "trait2"))
    expect_s3_class(plot_prior_posterior(other, prior = prior), "ggplot")
  }
})

test_that("LKJ overrides are validated and differences are explicit", {
  fit <- prior_test_fit()
  expect_warning(plot_prior_posterior(fit, lkj_prior = 1, n_prior = 50),
                 "differs")
  expect_warning(plot_forest(fit, lkj_prior = 1, n_prior = 50), "differs")
  expect_no_warning(plot_prior_posterior(fit, lkj_prior = 3, n_prior = 50))
  expect_error(plot_forest(fit, lkj_prior = 0, n_prior = 50), "`lkj_prior`")
  expect_warning(plot_forest(prior_test_fit(1), lkj_prior = 3, n_prior = 50),
                 "ignored")
  expect_warning(plot_forest(fit, lkj_prior = 3,
                             prior = animal_prior_draws(2, 50)), "ignored")
})

test_that("trait indices and recorded names determine labels", {
  fit <- prior_test_fit()
  fit <- fit[, c("h2[2]", "h2[1]", "corrG[1,2]")]
  attr(fit, "y_names") <- c("length", "width")
  d <- animal_parameter_draws(fit)
  expect_equal(unique(d$value[d$parameter == "length"]), 0.1)
  expect_equal(unique(d$value[d$parameter == "width"]), 0.2)
  expect_true("length ~ width" %in% d$parameter)
  bad <- fit[, c("h2[2]", "corrG[1,2]"), drop = FALSE]
  attr(bad, "y_names") <- NULL
  expect_error(animal_parameter_draws(bad), "h2.*indices")
  expect_error(animal_parameter_draws(fit, c("same", "same")), "trait_names")
})

test_that("explicit prior labels, values and source must match the posterior", {
  fit <- prior_test_fit()
  prior <- animal_prior_draws(2, 30, trait_names = c("trait1", "trait2"))
  expect_s3_class(plot_forest(fit, prior = prior), "ggplot")
  wrong <- prior
  wrong$parameter[wrong$parameter == "trait1"] <- "other"
  expect_error(plot_prior_posterior(fit, prior = wrong), "parameter")
  wrong <- prior[prior$parameter != "trait2", ]
  expect_error(plot_forest(fit, prior = wrong), "parameter")
  wrong <- prior
  wrong$type[1] <- "Genetic correlation"
  expect_error(plot_forest(fit, prior = wrong), "parameter|type")
  wrong <- prior
  wrong$value[1] <- NA_real_
  expect_error(plot_forest(fit, prior = wrong), "value")
  wrong <- prior
  wrong$value[1] <- 1.5
  expect_error(plot_forest(fit, prior = wrong), "value")
  wrong <- prior
  wrong$source[1] <- "posterior"
  expect_error(plot_forest(fit, prior = wrong), "source")
})

test_that("LKJ priors match analytical correlations for more than two traits", {
  set.seed(105)
  d <- animal_prior_draws(4, n = 4000, lkj_prior = 2)
  for (label in unique(d$parameter[d$type == "Genetic correlation"])) {
    rho <- d$value[d$parameter == label]
    expect_lt(abs(mean(rho)), 0.035)
    expect_lt(abs(var(rho) - 1 / 7), 0.02)
  }
})

test_that("multivariate h2 and diagonal residual fits retain their actual prior", {
  skip_if_no_cmdstan()
  sim <- simulate_animal(n_fam = 10, size = 4, G = diag(2), E = diag(2),
                         seed = 78)
  colnames(sim$Y) <- c("length", "width")
  fit <- suppressWarnings(lmm_animal(
    sim$Y, sim$X, sim$A, parametrization = "h2",
    residual_structure = "diagonal", beta_prior = c(2, 8), lkj_prior = 4,
    iter = 300, warmup = 150, chains = 1, cores = 1, seed = 1,
    refresh = 0, show_messages = FALSE))
  set.seed(2)
  d <- stanimal:::.prior_posterior(fit, NULL, NULL, 4000, NULL)
  h2 <- d$value[d$source == "Prior" & d$type == "Heritability"]
  rho <- d$value[d$source == "Prior" & d$type == "Genetic correlation"]
  expect_lt(abs(mean(h2) - 0.2), 0.02)
  expect_lt(abs(var(rho) - 1 / 9), 0.02)
  expect_setequal(unique(d$parameter[d$type == "Heritability"]), colnames(sim$Y))
})

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

test_that("a factor fit reconstructs its prior or accepts an explicit one", {
  skip_if_no_cmdstan()
  sim <- simulate_animal(n_fam = 10, size = 4, G = diag(c(1, 0.5)),
                         E = diag(2), seed = 74)
  fit <- suppressWarnings(lmm_animal(sim$Y, sim$X, sim$A,
                                     genetic_structure = "factor",
                                     n_factors = 1, beta_prior = c(2, 8),
                                     hs_global_scale = 0.7, hs_slab_scale = 1.5,
                                     column_decay = 0.6, iter = 300, warmup = 150,
                                     chains = 1, cores = 1, seed = 1,
                                     refresh = 0, show_messages = FALSE))
  set.seed(76)
  d <- stanimal:::.prior_posterior(fit, NULL, NULL, 4000, NULL)
  h2 <- d$value[d$source == "Prior" & d$type == "Heritability"]
  expect_lt(abs(mean(h2) - 0.2), 0.02)
  expect_equal(attr(fit, "hs_slab_scale"), 1.5)
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
