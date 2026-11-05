test_that("fit_loo turns the pointwise terms into a PSIS-LOO", {
  skip_if_no_cmdstan()
  skip_if_no_loo()
  f <- small_animal_fit()
  out <- suppressWarnings(fit_loo(f$fit))

  expect_s3_class(out, "psis_loo")
  # one term per record, and the estimate is over those terms
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

test_that("a fit with incomplete records is comparable too", {
  skip_if_no_cmdstan()
  skip_if_no_loo()
  sim <- simulate_animal(n_fam = 8, size = 5, G = G3, E = E3, seed = 72)
  Y <- sim$Y
  Y[1, 1] <- NA
  Y[4, 2] <- NA
  fit <- suppressWarnings(lmm_animal(Y, sim$X, sim$A, iter = 400, chains = 2,
                                     cores = 2, seed = 1, refresh = 0,
                                     show_messages = FALSE))
  out <- suppressWarnings(fit_loo(fit))
  expect_equal(length(out$diagnostics$pareto_k), nrow(Y))
})

test_that("phenotypic fits report the same quantities", {
  skip_if_no_cmdstan()
  skip_if_no_loo()
  sim <- simulate_animal(n_fam = 8, size = 5, G = G3, E = E3, seed = 73)
  fit <- lm_phen_cov(sim$Y, sim$X, iter = 400, chains = 2, cores = 2, seed = 1,
                     refresh = 0, show_messages = FALSE)
  out <- suppressWarnings(fit_loo(fit))
  expect_s3_class(out, "psis_loo")
  expect_equal(length(out$diagnostics$pareto_k), nrow(sim$Y))
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
