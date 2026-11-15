# correlation_prior = "size" / "shared_size": the size-factor prior on the
# correlations of the unstructured h2 model (animalModelGH2Size.stan).

# Four traits on one strong size axis, in G and in E, in families of full sibs.
size_sim <- function(seed, n_fam = 40, size = 5) {
  set.seed(seed)
  K <- 4
  n <- n_fam * size
  R <- matrix(0.7, K, K); diag(R) <- 1
  G <- R * tcrossprod(sqrt(c(1, 0.8, 1.2, 0.9)))
  E <- (matrix(0.4, K, K) + diag(0.6, K)) * tcrossprod(sqrt(c(1, 1.1, 0.9, 1)))
  A <- make_A(n_fam, size)
  X <- cbind(Int = 1, x = stats::rnorm(n, 50, 10))
  Y <- X %*% rbind(c(10, -5, 2, 4), c(0.05, 0.1, -0.02, 0)) +
    t(chol(A)) %*% rmvn(n, G) + rmvn(n, E)
  colnames(Y) <- c("a", "b", "c", "d")
  list(Y = Y, X = X, A = A, G = G, E = E)
}

size_fit <- function(sim, prior, iter = 800, warmup = 400, ...) {
  lmm_animal(sim$Y, sim$X, sim$A, correlation_prior = prior,
             iter = iter, warmup = warmup, chains = 2, cores = 2, seed = 1,
             refresh = 0, show_messages = FALSE, ...)
}

# draws of a two-trait h2 fit with recorded metadata, for the prior rebuild
size_prior_test_fit <- function() {
  vars <- c("h2[1]", "h2[2]", "corrG[1,2]")
  fit <- posterior::as_draws_matrix(
    matrix(rep(c(0.1, 0.2, 0.3), each = 20), 20, dimnames = list(NULL, vars)))
  attr(fit, "genetic_structure") <- "unstructured"
  attr(fit, "residual_structure") <- "unstructured"
  attr(fit, "parametrization") <- "h2"
  attr(fit, "custom_model") <- FALSE
  attr(fit, "lkj_prior") <- 2
  attr(fit, "beta_prior") <- c(2, 8)
  attr(fit, "y_names") <- c("trait1", "trait2")
  fit
}

# loadings on the correlation scale back to l = mu + sigma z
raw_loading <- function(s) s / sqrt(1 - s^2)

test_that("the size priors are refused where they are not defined", {
  sim <- simulate_animal(n_fam = 5, size = 4, G = G3, E = E3, seed = 3)
  call <- function(...) lmm_animal(sim$Y, sim$X, sim$A, iter = 10, ...)
  expect_error(lmm_animal(sim$Y[, 1], sim$X, sim$A, correlation_prior = "size"),
               "at least two traits")
  expect_error(call(correlation_prior = "size", genetic_structure = "factor",
                    n_factors = 1),
               "genetic_structure = \"unstructured\"")
  expect_error(call(correlation_prior = "shared_size",
                    residual_structure = "diagonal"),
               "residual_structure = \"unstructured\"")
  expect_error(call(correlation_prior = "size", model_file = "x.stan"),
               "cannot be combined with `model_file`")
  expect_error(call(correlation_prior = "size", parametrization = "variance"),
               "only available with `parametrization = \"h2\"`")
  expect_error(call(correlation_prior = "uniform"), "should be one of")
})

test_that("a size fit uses the h2 model, records its prior and reports the size axis", {
  skip_if_no_cmdstan()
  sim <- size_sim(seed = 11)
  warns <- character()
  fit <- withCallingHandlers(
    size_fit(sim, "size"),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  # lkj_prior and beta_prior are both used here, and parametrization was set
  expect_false(any(grepl("ignored", warns)))
  expect_clean_fit(fit)
  expect_equal(attr(fit, "correlation_prior"), "size")
  expect_equal(attr(fit, "parametrization"), "h2")

  vars <- fit$metadata()$stan_variables
  expect_true(all(c("rho_G", "rho_E", "sigma_size_G", "sigma_size_E",
                    "size_loadings_G", "size_loadings_E", "integration_G",
                    "integration_E", "z_size_G", "z_size_E") %in% vars))
  dm <- as.matrix(fit$draws(format = "draws_matrix"))
  expect_true(all(abs(dm[, grep("^size_loadings", colnames(dm))]) < 1))
  integ <- dm[, c("integration_G", "integration_E")]
  expect_true(all(integ >= 0 & integ <= 1))
  # h2 is still the diagonal ratio
  expect_equal(unname(dm[, sprintf("h2[%d]", 1:4)]),
               unname(dm[, sprintf("G[%d,%d]", 1:4, 1:4)] /
                        dm[, sprintf("P[%d,%d]", 1:4, 1:4)]),
               tolerance = 1e-6)
  # the genetic correlations are all 0.7, so the size axis must be found
  expect_gt(mean(dm[, "rho_G"]), 0.3)

  # the draw tools treat it like any other unstructured fit
  d <- fit_draws(fit)
  expect_equal(dim(d$size_loadings_G), c(4L, nrow(dm)))
  expect_equal(rownames(d$size_loadings_G), colnames(sim$Y))
  expect_equal(dim(breeding_values(fit, se = FALSE)$values), dim(sim$Y))
  skip_if_not_installed("loo")
  expect_s3_class(suppressWarnings(fit_loo(fit)), "psis_loo")
})

test_that("the QR construction is the correlation matrix D (l l' + R) D", {
  skip_if_no_cmdstan()
  sim <- size_sim(seed = 12, n_fam = 20)
  fit <- suppressWarnings(size_fit(sim, "size", iter = 200, warmup = 100))
  dm <- as.matrix(fit$draws(c("rho_G", "sigma_size_G", "z_size_G", "L_R_G",
                              "corrG"), format = "draws_matrix"))
  dm <- unclass(dm)
  K <- 4
  err <- vapply(seq_len(20), function(i) {
    l <- sqrt(dm[i, "rho_G"] / (1 - dm[i, "rho_G"])) +
      dm[i, "sigma_size_G"] * dm[i, sprintf("z_size_G[%d]", 1:K)]
    L <- matrix(dm[i, sprintf("L_R_G[%d,%d]", rep(1:K, K), rep(1:K, each = K))],
                K, K)
    d <- 1 / sqrt(1 + l^2)
    Omega <- (tcrossprod(l) + tcrossprod(L)) * tcrossprod(d)
    corrG <- matrix(dm[i, sprintf("corrG[%d,%d]", rep(1:K, K), rep(1:K, each = K))],
                    K, K)
    max(abs(Omega - corrG))
  }, numeric(1))
  expect_lt(max(err), 1e-5)
})

test_that("shared_size gives G and E one standardized loading pattern", {
  skip_if_no_cmdstan()
  sim <- size_sim(seed = 13)
  fit <- suppressWarnings(size_fit(sim, "shared_size"))
  expect_clean_fit(fit)
  expect_equal(attr(fit, "correlation_prior"), "shared_size")
  vars <- fit$metadata()$stan_variables
  expect_true("z_size_G" %in% vars)
  dm <- unclass(as.matrix(fit$draws(format = "draws_matrix")))
  expect_false(any(grepl("^z_size_E", colnames(dm))))
  # (l_G - mu_G) / sigma_G and (l_E - mu_E) / sigma_E are the same vector z
  z <- function(m) {
    s <- dm[, sprintf("size_loadings_%s[%d]", m, 1:4)]
    rho <- dm[, if (m == "G") "rho_G" else "rho_E[1]"]
    sigma <- dm[, if (m == "G") "sigma_size_G" else "sigma_size_E[1]"]
    (raw_loading(s) - sqrt(rho / (1 - rho))) / sigma
  }
  keep <- dm[, "sigma_size_G"] > 0.05 & dm[, "sigma_size_E[1]"] > 0.05
  expect_lt(max(abs(z("G")[keep, ] - z("E")[keep, ])), 1e-3)
})

test_that("the size prior draws are the model's prior", {
  set.seed(5)
  p <- animal_prior_draws(4, 4000, lkj_prior = 2, beta_prior = c(2, 8),
                          correlation_prior = "size")
  h2 <- p$value[p$type == "Heritability"]
  r <- p$value[p$type == "Genetic correlation"]
  expect_lt(abs(mean(h2) - 0.2), 0.01)
  # the same transformation as the Stan file, by direct simulation: rho, sigma,
  # z and an LKJ(2) remainder
  set.seed(6)
  ref <- replicate(4000, {
    rho <- rbeta(1, 1, 1)
    l <- sqrt(rho / (1 - rho)) + abs(rnorm(1)) * rnorm(4)
    R <- stanimal:::.rlkj(4, 2)
    d <- 1 / sqrt(1 + l^2)
    ((tcrossprod(l) + R) * tcrossprod(d))[1, 2]
  })
  expect_lt(abs(mean(r) - mean(ref)), 0.03)
  expect_lt(abs(sd(r) - sd(ref)), 0.03)
  # and it is not the LKJ prior, which is centred on zero
  expect_gt(mean(r), 0.3)
  # shared_size has the same marginal prior on G
  set.seed(5)
  q <- animal_prior_draws(4, 4000, lkj_prior = 2, beta_prior = c(2, 8),
                          correlation_prior = "shared_size")
  expect_equal(p$value, q$value)
  # the size priors only exist with the h2 parametrization: Beta(1, 1) default
  set.seed(7)
  flat <- animal_prior_draws(3, 4000, correlation_prior = "size")
  expect_lt(abs(mean(flat$value[flat$type == "Heritability"]) - 0.5), 0.02)
  expect_error(animal_prior_draws(1, correlation_prior = "size"),
               "at least two traits")
})

test_that("plot_prior_posterior rebuilds the size prior a fit recorded", {
  fit <- size_prior_test_fit()
  attr(fit, "correlation_prior") <- "size"
  set.seed(8)
  d <- stanimal:::.prior_posterior(fit, NULL, NULL, 4000, NULL)
  r <- d$value[d$source == "Prior" & d$type == "Genetic correlation"]
  expect_gt(mean(r), 0.3)
  # without the attribute (an older fit) it is the LKJ prior
  attr(fit, "correlation_prior") <- NULL
  d <- stanimal:::.prior_posterior(fit, NULL, NULL, 4000, NULL)
  r <- d$value[d$source == "Prior" & d$type == "Genetic correlation"]
  expect_lt(abs(mean(r)), 0.03)
  # an unrecognized value is not guessed
  attr(fit, "correlation_prior") <- "uniform"
  expect_error(stanimal:::.prior_posterior(fit, NULL, NULL, 100, NULL),
               "cannot be reconstructed")
})

test_that("plot_prior_posterior works on a sampled size fit", {
  skip_if_no_cmdstan()
  sim <- size_sim(seed = 14, n_fam = 20)
  fit <- suppressWarnings(size_fit(sim, "shared_size", iter = 200, warmup = 100,
                                   beta_prior = c(2, 8)))
  set.seed(9)
  d <- stanimal:::.prior_posterior(fit, NULL, NULL, 2000, NULL)
  prior <- d[d$source == "Prior", ]
  expect_lt(abs(mean(prior$value[prior$type == "Heritability"]) - 0.2), 0.02)
  expect_gt(mean(prior$value[prior$type == "Genetic correlation"]), 0.3)
  expect_s3_class(plot_prior_posterior(fit, n_prior = 200), "ggplot")
})

test_that("the Laplace approximation takes the size priors", {
  skip_if_no_cmdstan()
  sim <- size_sim(seed = 15, n_fam = 20)
  fit <- suppressWarnings(lmm_animal_laplace(sim$Y, sim$X, sim$A,
                                             correlation_prior = "shared_size",
                                             draws = 200, seed = 1))
  q <- fit$summary$quantity
  expect_true(all(c("rho_G", "rho_E[1]", "integration_G", "size_loadings_G[a]",
                    "size_loadings_E[d]") %in% q))
  expect_error(lmm_animal_laplace(sim$Y, sim$X, sim$A,
                                  correlation_prior = "size",
                                  parametrization = "variance"),
               "only available")
})
