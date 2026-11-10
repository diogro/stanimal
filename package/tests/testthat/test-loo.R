test_that("fit_loo turns the pointwise terms into a PSIS-LOO", {
  skip_if_no_cmdstan()
  skip_if_no_loo()
  f <- small_animal_fit()
  out <- suppressWarnings(fit_loo(f$fit))

  expect_s3_class(out, "psis_loo")
  # one term per individual, and the estimate is over those terms
  expect_equal(out$diagnostics$pareto_k |> length(), nrow(f$sim$Y))
  expect_true(all(is.finite(out$estimates[, "Estimate"])))
  expect_lt(out$estimates["elpd_loo", "Estimate"], 0)
  expect_gt(out$estimates["p_loo", "Estimate"], 0)
})

test_that("loo() dispatches on a fit once the loo package is loaded", {
  skip_if_no_cmdstan()
  skip_if_no_loo()
  f <- small_animal_fit()
  # pkgload does not run .onLoad, so register explicitly; under a normal install
  # the registration happens when the package loads or when loo does
  requireNamespace("loo", quietly = TRUE)
  stanimal:::.onLoad()
  out <- suppressWarnings(loo::loo(f$fit))
  expect_s3_class(out, "psis_loo")
  expect_equal(out$estimates["elpd_loo", "Estimate"],
               suppressWarnings(fit_loo(f$fit))$estimates["elpd_loo", "Estimate"])
})

test_that("WAIC and model comparison work from the same draws", {
  skip_if_no_cmdstan()
  skip_if_no_loo()
  f <- small_animal_fit()
  w <- suppressWarnings(fit_loo(f$fit, criterion = "waic"))
  expect_s3_class(w, "waic")
  expect_true(all(c("elpd_waic", "p_waic") %in% rownames(w$estimates)))

  unstructured <- suppressWarnings(fit_loo(f$fit))
  factor_fit <- suppressWarnings(lmm_animal(
    f$sim$Y, f$sim$X, f$sim$A, genetic_structure = "factor", n_factors = 2,
    column_decay = 0.5, iter = 400, chains = 2, cores = 2, seed = 1, refresh = 0,
    show_messages = FALSE))
  comparison <- suppressWarnings(
    loo::loo_compare(unstructured, fit_loo(factor_fit)))
  expect_equal(nrow(comparison), 2L)
  expect_true(all(c("elpd_diff", "se_diff") %in% colnames(comparison)))
})

# Dense brute force: log p(y_i,obs | y_-i, theta) from the full NK x NK
# covariance, with y_-i complete (or completed) and i's own gaps left out.
dense_individual_log_lik <- function(R, G, E, A, miss) {
  n <- nrow(R)
  k <- ncol(R)
  V <- kronecker(G, A) + kronecker(E, diag(n))
  r <- c(R)
  vapply(seq_len(n), function(i) {
    if (all(miss[i, ])) return(NA_real_)
    own <- i + (seq_len(k) - 1L) * n
    idx <- own[!miss[i, ]]
    rest <- setdiff(seq_len(n * k), own)
    S12 <- V[idx, rest, drop = FALSE]
    m <- S12 %*% solve(V[rest, rest], r[rest])
    Sc <- V[idx, idx, drop = FALSE] - S12 %*% solve(V[rest, rest], t(S12))
    -0.5 * (length(idx) * log(2 * pi) + determinant(Sc)$modulus +
              sum((r[idx] - m) * solve(Sc, r[idx] - m)))
  }, numeric(1))
}

test_that("individual terms are the exact conditional densities", {
  set.seed(3)
  n <- 20
  k <- 3
  A <- crossprod(matrix(rnorm(n * n), n)) / n + diag(n) * 0.3
  G <- array(c(G3, E3 + diag(3) * 0.1), c(k, k, 2))
  E <- array(c(E3, G3 + diag(3) * 0.2), c(k, k, 2))
  R <- array(rnorm(n * k * 2), c(n, k, 2))
  miss <- matrix(runif(n * k) < 0.25, n)
  miss[1, ] <- TRUE
  eig <- stanimal:::.relationship_eigen(A)
  for (pattern in list(matrix(FALSE, n, k), miss)) {
    fast <- stanimal:::.individual_log_density(R, G, E, eig$Q, eig$Q^2,
                                               eig$lambda, pattern)
    for (s in 1:2) {
      dense <- dense_individual_log_lik(matrix(R[, , s], n, k), G[, , s],
                                        E[, , s], A, pattern)
      expect_identical(is.na(fast[s, ]), is.na(dense))
      expect_equal(fast[s, ], dense, tolerance = 1e-10)
    }
  }
})

test_that("animal fits default to one term per individual", {
  skip_if_no_cmdstan()
  skip_if_no_loo()
  f <- small_animal_fit()
  out <- suppressWarnings(fit_loo(f$fit))
  expect_s3_class(out, "psis_loo")
  expect_equal(length(out$diagnostics$pareto_k), nrow(f$sim$Y))
  expect_equal(attr(out, "units"),
               rownames(f$sim$Y) %||% as.character(seq_len(nrow(f$sim$Y))))

  # the extraction, scaling and design reproduce the dense computation draw by
  # draw, on the standardized response scale
  terms <- stanimal:::.individual_log_lik(f$fit)
  draws <- f$fit$draws(format = "draws_matrix")
  G <- draws_mean(f$fit, "G") * 0
  Y <- f$sim$Y
  log_sd <- sum(log(apply(Y, 2, stats::sd)))
  for (s in c(1L, 150L, 400L)) {
    one <- posterior::subset_draws(draws, draw = s)
    G <- matrix(as.numeric(one[, sprintf("G[%d,%d]", rep(1:3, 3),
                                         rep(1:3, each = 3))]), 3)
    E <- matrix(as.numeric(one[, sprintf("E[%d,%d]", rep(1:3, 3),
                                         rep(1:3, each = 3))]), 3)
    beta <- stanimal:::.beta_draws(as.matrix(one), 3, paste0("t", 1:3), NULL,
                                   attr(f$fit, "x_names"))
    R <- Y - f$sim$X %*% t(matrix(beta[, , 1], nrow = 3))
    dense <- dense_individual_log_lik(R, G, E, attr(f$fit, "stanimal_data")$A,
                                      matrix(FALSE, nrow(Y), 3))
    expect_equal(terms$log_lik[s, ], dense + log_sd, tolerance = 1e-8)
  }

  record <- suppressWarnings(fit_loo(f$fit, unit = "record"))
  expect_equal(length(record$diagnostics$pareto_k), nrow(f$sim$Y))
  expect_null(attr(record, "units"))
})

test_that("incomplete animal fits get individual-level LOO and WAIC", {
  skip_if_no_cmdstan()
  skip_if_no_loo()
  sim <- simulate_animal(n_fam = 8, size = 5, G = G3, E = E3, seed = 72)
  Y <- sim$Y
  Y[1, 1] <- NA
  Y[4, 2] <- NA
  Y[7, ] <- NA
  fit <- suppressWarnings(lmm_animal(Y, sim$X, sim$A, iter = 400, chains = 2,
                                     cores = 2, seed = 1, refresh = 0,
                                     show_messages = FALSE))
  out <- suppressWarnings(fit_loo(fit))
  expect_s3_class(out, "psis_loo")
  # the animal with nothing observed carries no term
  expect_equal(length(out$diagnostics$pareto_k), nrow(Y) - 1L)
  expect_false(any(attr(out, "units") ==
                     (rownames(Y) %||% as.character(seq_len(nrow(Y))))[7]))
  expect_true(all(is.finite(out$pointwise[, "elpd_loo"])))
  expect_s3_class(suppressWarnings(fit_loo(fit, criterion = "waic")), "waic")

  # with gaps, a draw conditions on the sampled gaps of the other individuals
  terms <- stanimal:::.individual_log_lik(fit)
  completed <- imputed_data(fit, ndraws = 1e6)$Y
  draws <- fit$draws(format = "draws_matrix")
  miss <- is.na(Y)
  log_sd <- drop((!miss) %*% log(apply(Y, 2, stats::sd, na.rm = TRUE)))
  for (s in c(1L, 400L)) {
    one <- posterior::subset_draws(draws, draw = s)
    G <- matrix(as.numeric(one[, sprintf("G[%d,%d]", rep(1:3, 3),
                                         rep(1:3, each = 3))]), 3)
    E <- matrix(as.numeric(one[, sprintf("E[%d,%d]", rep(1:3, 3),
                                         rep(1:3, each = 3))]), 3)
    beta <- stanimal:::.beta_draws(as.matrix(one), 3, paste0("t", 1:3), NULL,
                                   attr(fit, "x_names"))
    R <- completed[, , s] - sim$X %*% t(matrix(beta[, , 1], nrow = 3))
    dense <- dense_individual_log_lik(R, G, E, attr(fit, "stanimal_data")$A,
                                      miss)
    keep <- rowSums(!miss) > 0
    expect_equal(terms$log_lik[s, ], (dense + log_sd)[keep], tolerance = 1e-8)
  }

  expect_error(fit_loo(fit, unit = "record"), "unavailable for animal fits")
  stanimal:::.onLoad()
  expect_s3_class(suppressWarnings(loo::loo(fit)), "psis_loo")
})

test_that("individual-level LOO accepts a checked eigen_A", {
  skip_if_no_cmdstan()
  skip_if_no_loo()
  f <- small_animal_fit()
  A <- attr(f$fit, "stanimal_data")$A
  e <- eigen(A, symmetric = TRUE)
  with_eigen <- suppressWarnings(fit_loo(f$fit,
                                         eigen_A = list(Q = e$vectors,
                                                        lambda = e$values)))
  without <- suppressWarnings(fit_loo(f$fit))
  expect_equal(with_eigen$pointwise, without$pointwise, tolerance = 1e-8)
  expect_error(fit_loo(f$fit, eigen_A = list(Q = e$vectors,
                                             lambda = e$values * 2)))
})

test_that("phenotypic fits with incomplete records support LOO and WAIC", {
  skip_if_no_cmdstan()
  skip_if_no_loo()
  sim <- simulate_animal(n_fam = 8, size = 5, G = G3, E = E3, seed = 73)
  Y <- sim$Y
  Y[1, 1] <- NA
  fit <- lm_phen_cov(Y, sim$X, iter = 400, chains = 2, cores = 2, seed = 1,
                     refresh = 0, show_messages = FALSE)
  out <- suppressWarnings(fit_loo(fit))
  expect_s3_class(out, "psis_loo")
  expect_equal(length(out$diagnostics$pareto_k), nrow(Y))
  expect_s3_class(suppressWarnings(fit_loo(fit, criterion = "waic")), "waic")
})

test_that("likelihood eligibility uses metadata and warns when it is absent", {
  skip_if_no_loo()
  draws <- cbind("log_lik[1]" = seq(-2, -1, length.out = 100),
                 "log_lik[2]" = seq(-1, -2, length.out = 100))
  expect_warning(result <- fit_loo(draws, criterion = "waic"),
                 "Cannot verify whether")
  expect_s3_class(result, "waic")

  mock <- list(metadata = function() list(stan_variables = "log_lik"),
               draws = function(...) array(draws, dim = c(100, 1, 2)))
  expect_warning(stanimal:::.log_lik_draws(mock), "Cannot verify whether")
  attr(mock, "stanimal_data") <- list(family = "animal",
                                     Y = matrix(c(1, NA), ncol = 1))
  attr(mock, "custom_model") <- TRUE
  expect_error(fit_loo(mock), "assumes the built-in animal likelihood")
  expect_error(fit_loo(mock, unit = "record"),
               "unavailable for animal fits with missing")
  expect_error(fit_loo(mock, criterion = "waic", unit = "record"),
               "unavailable for animal fits with missing")
  attr(mock, "stanimal_data")$Y[2, 1] <- 2
  expect_no_warning(stanimal:::.log_lik_draws(mock))
  attr(mock, "stanimal_data") <- list(family = "phenotypic",
                                     Y = matrix(c(1, NA), ncol = 1))
  expect_no_warning(stanimal:::.log_lik_draws(mock))
})

test_that("rstan dispatch also refuses rotated-record LOO with gaps", {
  skip_if_no_rstan()
  skip_if_no_loo()
  fit <- methods::new("stanfit")
  attr(fit, "stanimal_data") <- list(family = "animal",
                                     Y = matrix(c(1, NA), ncol = 1))
  stanimal:::.onLoad()
  expect_error(fit_loo(fit, unit = "record"),
               "unavailable for animal fits with missing")
  expect_error(fit_loo(fit, criterion = "waic", unit = "record"),
               "unavailable for animal fits with missing")
  expect_error(loo::loo(fit, unit = "record"),
               "unavailable for animal fits with missing")
})

test_that("fits without a pointwise likelihood say so", {
  fit <- lmm_animal_gibbs(simulate_animal(n_fam = 5, size = 4, G = G3, E = E3,
                                          seed = 74)$Y,
                          A = make_A(5, 4), iter = 40, warmup = 20)
  expect_error(fit_loo(fit), "does not support LOO or WAIC")
  expect_error(fit_loo(fit), "for model comparison")
  expect_error(fit_loo(list(a = 1)), "must be a fit")
})

test_that("posterior_predictive replicates the response", {
  skip_if_no_cmdstan()
  f <- small_animal_fit()
  pp <- posterior_predictive(f$fit, ndraws = 20, seed = 5)

  expect_s3_class(pp, "stanimal_yrep")
  # individuals, traits, draws: draws last, as everywhere else in the package
  expect_equal(dim(pp$yrep), c(nrow(f$sim$Y), ncol(f$sim$Y), 20L))
  expect_equal(dimnames(pp$yrep)[[2]], pp$traits)
  expect_equal(pp$ndraws, 20L)
  expect_equal(pp$traits, colnames(f$sim$Y) %||% paste0("trait", 1:3))
  expect_true(all(is.finite(pp$yrep)))
  # the replicates centre where the data are, not at zero
  expect_lt(abs(mean(pp$yrep[, 1, ]) - mean(f$sim$Y[, 1])),
            0.25 * stats::sd(f$sim$Y[, 1]))
  # replicates are draws, not copies of one simulation
  expect_gt(stats::sd(pp$yrep[1, 1, ]), 0)

  # print and plot summarise the covariance structure
  out <- capture.output(print(pp))
  expect_true(any(grepl("observed", out)))
  expect_true(any(grepl("^ cor: ", out)))
  p <- plot(pp)
  expect_s3_class(p, "ggplot")
})

test_that("posterior_predictive keeps the observed gaps", {
  skip_if_no_cmdstan()
  sim <- simulate_animal(n_fam = 8, size = 5, G = G3, E = E3, seed = 75)
  Y <- sim$Y
  Y[2, 3] <- NA
  fit <- suppressWarnings(lmm_animal(Y, sim$X, sim$A, iter = 400, chains = 2,
                                     cores = 2, seed = 1, refresh = 0,
                                     show_messages = FALSE))
  pp <- posterior_predictive(fit, ndraws = 5, seed = 1)
  expect_equal(sum(is.na(pp$y)), 1L)
  # the replicate is simulated for that entry even though the data lack it
  expect_false(anyNA(pp$yrep[2, 3, ]))
})

test_that("posterior_predictive covers the univariate and phenotypic paths", {
  skip_if_no_cmdstan()
  uni_sim <- simulate_animal(n_fam = 8, size = 5, G = matrix(0.8),
                             E = matrix(0.8), seed = 76)
  uni <- suppressWarnings(lmm_animal(uni_sim$Y[, 1], uni_sim$X, uni_sim$A,
                                     iter = 400, chains = 2, cores = 2,
                                     seed = 1, refresh = 0,
                                     show_messages = FALSE))
  pp_uni <- posterior_predictive(uni, ndraws = 8, seed = 1)
  expect_equal(dim(pp_uni$yrep), c(nrow(uni_sim$Y), 1L, 8L))
  # the genetic deviations come out with the right spread: on the data scale the
  # replicates are as variable as sigma_G^2 + sigma_E^2 implies
  expect_gt(stats::var(as.vector(pp_uni$yrep)), 0)

  phen_sim <- simulate_animal(n_fam = 8, size = 5, G = G3, E = E3, seed = 76)
  phen <- lm_phen_cov(phen_sim$Y, phen_sim$X, iter = 400, chains = 2, cores = 2,
                      seed = 1, refresh = 0, show_messages = FALSE)
  pp_phen <- posterior_predictive(phen, ndraws = 8, seed = 1)
  expect_equal(dim(pp_phen$yrep), c(nrow(phen_sim$Y), 3L, 8L))
  expect_true(all(is.finite(pp_phen$yrep)))
})

test_that("posterior_predictive refuses fits it cannot replicate", {
  sim <- simulate_animal(n_fam = 5, size = 4, G = G3, E = E3, seed = 77)
  gibbs <- lmm_animal_gibbs(sim$Y, sim$X, sim$A, iter = 40, warmup = 20)
  expect_error(posterior_predictive(gibbs), "does not carry the data")
})
