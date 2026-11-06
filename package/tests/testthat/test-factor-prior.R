factor_prior_fit <- function(parametrization = "h2", residual = "unstructured") {
  vars <- c("h2[1]", "h2[2]", "h2[3]", "corrG[1,2]",
            "corrG[1,3]", "corrG[2,3]")
  fit <- posterior::as_draws_matrix(matrix(rep(0.2, 20 * length(vars)), 20,
                                           dimnames = list(NULL, vars)))
  attr(fit, "genetic_structure") <- "factor"
  attr(fit, "residual_structure") <- residual
  attr(fit, "parametrization") <- parametrization
  attr(fit, "custom_model") <- FALSE
  attr(fit, "n_factors") <- 2L
  attr(fit, "hs_global_scale") <- 0.7
  attr(fit, "hs_slab_scale") <- 1.5
  attr(fit, "column_decay") <- 0.6
  attr(fit, "beta_prior") <- c(2, 8)
  attr(fit, "y_names") <- c("a", "b", "c")
  fit
}

test_that("factor G priors are reconstructed for both residual structures", {
  for (residual in c("unstructured", "diagonal")) {
    fit <- factor_prior_fit(residual = residual)
    set.seed(81)
    d <- stanimal:::.prior_posterior(fit, NULL, NULL, 4000, NULL)
    set.seed(81)
    expect_equal(d, stanimal:::.prior_posterior(fit, NULL, NULL, 4000, NULL))
    prior <- d[d$source == "Prior", ]
    expect_setequal(unique(as.character(prior$parameter)),
                    c("a", "b", "c", "a ~ b", "a ~ c", "b ~ c"))
    h2 <- prior$value[prior$type == "Heritability"]
    expect_lt(abs(mean(h2) - 0.2), 0.02)
    expect_true(all(h2 >= 0 & h2 <= 1))
    expect_true(all(abs(prior$value[prior$type == "Genetic correlation"]) <= 1))
    expect_s3_class(plot_prior_posterior(fit, n_prior = 100), "ggplot")
    expect_s3_class(plot_forest(fit, n_prior = 100), "ggplot")
    expect_warning(plot_prior_posterior(fit, lkj_prior = 2, n_prior = 100),
                   "ignored.*factor")
  }
})

test_that("factor priors refuse incomplete and unsupported configurations", {
  fit <- factor_prior_fit()
  for (field in c("n_factors", "hs_global_scale", "hs_slab_scale",
                  "column_decay", "beta_prior", "parametrization",
                  "residual_structure", "custom_model")) {
    incomplete <- fit
    attr(incomplete, field) <- NULL
    expect_error(plot_prior_posterior(incomplete, n_prior = 30),
                 "cannot.*reconstruct")
  }
  for (field in c("n_factors", "hs_global_scale", "hs_slab_scale",
                  "column_decay", "beta_prior")) {
    bad <- fit
    attr(bad, field) <- 0
    expect_error(plot_forest(bad, n_prior = 30), field)
  }
  bad <- fit
  attr(bad, "n_factors") <- 3L
  expect_error(plot_forest(bad, n_prior = 30), "n_factors")
  attr(fit, "residual_structure") <- "shared_factor"
  expect_error(plot_forest(fit, n_prior = 30), "shared.factor")
  attr(fit, "custom_model") <- TRUE
  expect_error(plot_forest(fit, n_prior = 30), "custom")
})

test_that("the historical factor variance parametrization uses its own h2 prior", {
  fit <- factor_prior_fit("variance")
  attr(fit, "beta_prior") <- NULL
  set.seed(82)
  d <- stanimal:::.prior_posterior(fit, NULL, NULL, 4000, NULL)
  h2 <- d$value[d$source == "Prior" & d$type == "Heritability"]
  expect_lt(abs(mean(h2) - 0.5), 0.02)
  expect_lt(abs(var(h2) - 0.125), 0.01)
})

test_that("factor count, horseshoe scales and decay affect the generated prior", {
  fit <- factor_prior_fit()
  set.seed(84)
  reference <- stanimal:::.prior_posterior(fit, NULL, NULL, 4000, NULL)
  rho <- reference$value[reference$source == "Prior" &
                           reference$parameter == "a ~ b"]
  for (field in c("n_factors", "hs_global_scale", "hs_slab_scale",
                  "column_decay")) {
    other <- fit
    attr(other, field) <- switch(field, n_factors = 1L,
                                 hs_global_scale = 0.05, hs_slab_scale = 0.2,
                                 column_decay = 0.05)
    set.seed(84)
    d <- stanimal:::.prior_posterior(other, NULL, NULL, 4000, NULL)
    alternative <- d$value[d$source == "Prior" & d$parameter == "a ~ b"]
    expect_gt(abs(var(rho) - var(alternative)), 0.002)
  }
})

test_that("simulated factor correlations are positive definite", {
  set.seed(83)
  for (i in seq_len(200)) {
    R <- stanimal:::.factor_g_prior_correlation(4, 3, 0.7, 1.5, 0.6)
    expect_equal(diag(R), rep(1, 4), tolerance = 1e-12)
    expect_equal(R, t(R), tolerance = 1e-12)
    expect_gt(min(eigen(R, symmetric = TRUE, only.values = TRUE)$values), 0)
  }
  for (slab in c(1e-200, 1e200)) {
    R <- stanimal:::.factor_g_prior_correlation(4, 3, 0.7, slab, 0.6)
    expect_true(all(is.finite(R)))
    expect_equal(diag(R), rep(1, 4), tolerance = 1e-12)
  }
})

test_that("factor G prior simulation agrees with the Stan transformation", {
  skip_if_no_cmdstan()
  stan <- readLines(system.file("stan", "animalModelGFactor.stan",
                                package = "stanimal"))
  start <- grep("  // Regularized horseshoe \\(Piironen", stan, fixed = FALSE)
  end <- grep("^model \\{$", stan)
  expect_length(start, 1)
  expect_length(end, 1)
  transform <- stan[seq.int(start, end - 2L)]
  code <- paste(c(
    "data { int<lower=2> K; int<lower=1> M; real<lower=0> hs_global_scale;",
    "real<lower=0> hs_slab_scale; real<lower=0,upper=1> column_decay;",
    "vector<lower=0>[2] beta_prior; }",
    "generated quantities {",
    "matrix[K,M] loadings_raw; matrix[K,M] local_scale;",
    "real global_scale = abs(cauchy_rng(0, hs_global_scale));",
    "vector[M] column_scale; matrix[K,M] loadings_corr;",
    "vector[K] uniqueness_corr; vector[K] communality_corr;",
    "vector[K] h2; matrix[K,K] corrG;",
    "for (m in 1:M) column_scale[m] = pow(column_decay, m-1);",
    "for (k in 1:K) { h2[k] = beta_rng(beta_prior[1], beta_prior[2]);",
    "for (m in 1:M) { loadings_raw[k,m] = normal_rng(0,1);",
    "local_scale[k,m] = abs(cauchy_rng(0,1)); } }",
    transform,
    "corrG = loadings_corr * loadings_corr' + diag_matrix(uniqueness_corr);",
    "}"), collapse = "\n")
  model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code),
                                  quiet = TRUE)
  settings <- list(c(1, 0.4, 0.7, 0.5), c(2, 2, 2, 0.9))
  for (s in settings) {
    data <- list(K = 3L, M = as.integer(s[1]), hs_global_scale = s[2],
                 hs_slab_scale = s[3], column_decay = s[4],
                 beta_prior = c(2, 8))
    reference <- model$sample(data = data, fixed_param = TRUE, chains = 1,
                              iter_warmup = 0, iter_sampling = 5000,
                              seed = 91, refresh = 0, show_messages = FALSE)
    ref <- as.matrix(reference$draws(c("h2", "corrG"), format = "draws_matrix"))
    fit <- factor_prior_fit()
    attr(fit, "n_factors") <- data$M
    attr(fit, "hs_global_scale") <- data$hs_global_scale
    attr(fit, "hs_slab_scale") <- data$hs_slab_scale
    attr(fit, "column_decay") <- data$column_decay
    set.seed(92)
    d <- stanimal:::.prior_posterior(fit, NULL, NULL, 5000, NULL)
    prior <- d[d$source == "Prior", ]
    for (pair in list(c("a", "h2[1]"), c("a ~ b", "corrG[1,2]"),
                      c("a ~ c", "corrG[1,3]"))) {
      r <- prior$value[prior$parameter == pair[1]]
      expected <- ref[, pair[2]]
      expect_lt(max(abs(quantile(r, c(0.1, 0.5, 0.9)) -
                        quantile(expected, c(0.1, 0.5, 0.9)))), 0.06)
      expect_lt(abs(mean(r) - mean(expected)), 0.04)
      expect_lt(unname(stats::ks.test(r, expected, exact = FALSE)$statistic), 0.04)
    }
  }
})
