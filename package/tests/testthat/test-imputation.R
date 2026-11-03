# Completed datasets: the response with its gaps filled, one per draw.

test_that("a fit without gaps has nothing to impute", {
  sim <- simulate_animal(n_fam = 5, size = 4, G = G3, E = E3, seed = 91)
  fit <- lm_phen_cov(sim$Y, sim$X, iter = 300, chains = 2, cores = 2, seed = 1,
                     refresh = 0, show_messages = FALSE)
  expect_error(imputed_data(fit), "no missing values")
})

test_that("an animal fit hands back the gaps it sampled", {
  skip_if_no_cmdstan()
  sim <- simulate_animal(n_fam = 20, size = 8, G = G3, E = E3, seed = 92)
  truth <- sim$Y
  Y <- truth
  Y[3, 1] <- NA
  Y[7, 2] <- NA
  Y[40, 2] <- NA
  fit <- suppressWarnings(lmm_animal(Y, sim$X, sim$A, iter = 800, chains = 2,
                                     cores = 2, seed = 1, refresh = 0,
                                     show_messages = FALSE))
  imp <- imputed_data(fit, ndraws = 200, seed = 1)

  expect_s3_class(imp, "stanimal_imputation")
  # individuals by traits by draws, as everywhere else
  expect_equal(dim(imp$Y), c(nrow(Y), ncol(Y), 200L))
  expect_equal(dimnames(imp$Y)[[2]], imp$traits)
  expect_equal(imp$gaps$individual, c(3L, 7L, 40L))
  expect_equal(imp$gaps$trait, c(1L, 2L, 2L))

  # the observed cells are the data in every draw, and nothing is left missing
  expect_false(anyNA(imp$Y))
  expect_equal(imp$Y[, , 1][!is.na(Y)], Y[!is.na(Y)])
  expect_equal(imp$Y[, , 200][!is.na(Y)], Y[!is.na(Y)])
  # the gaps do vary: these are draws, not one fill repeated
  gap <- imp$Y[3, 1, ]
  expect_gt(stats::sd(gap), 0)
  # and they are on the data scale, not the standardized one the parameters are
  # declared on: a standardized value would sit outside the observed range of a
  # trait whose mean is 13
  positions <- cbind(imp$gaps$individual, imp$gaps$trait)
  for (i in seq_len(nrow(positions))) {
    j <- positions[i, 2L]
    observed <- Y[, j]
    expect_true(all(imp$Y[positions[i, 1L], j, ] > min(observed, na.rm = TRUE) -
                      stats::sd(observed, na.rm = TRUE) &
                      imp$Y[positions[i, 1L], j, ] < max(observed, na.rm = TRUE) +
                      stats::sd(observed, na.rm = TRUE)))
  }
  # and they land near the value that was removed, within the spread of the trait
  means <- apply(imp$Y, c(1L, 2L), mean)[positions]
  error <- abs(means - truth[positions])
  expect_true(all(error < 1.5 * apply(truth, 2L, stats::sd)[positions[, 2L]]))

  out <- capture.output(print(imp))
  expect_true(any(grepl("completed datasets", out)))
  # printed in individual order, not the order the parameters were declared in
  expect_true(any(grepl("^ +3 ", out)))
})

test_that("a phenotypic fit draws each gap from its conditional distribution", {
  skip_if_no_cmdstan()
  sim <- simulate_animal(n_fam = 30, size = 8, G = G3, E = E3, seed = 93)
  Y <- sim$Y
  Y[5, 2] <- NA                    # traits 1 and 3 stay observed for that row
  fit <- suppressWarnings(lm_phen_cov(Y, sim$X, iter = 800, chains = 2,
                                      cores = 2, seed = 1, refresh = 0,
                                      show_messages = FALSE))
  imp <- imputed_data(fit, ndraws = 300, seed = 1)

  expect_equal(dim(imp$Y), c(nrow(Y), ncol(Y), 300L))
  expect_false(anyNA(imp$Y))
  expect_equal(imp$Y[, , 1][!is.na(Y)], Y[!is.na(Y)])

  # conditioning on two correlated traits must beat the marginal spread: if the
  # observed traits were ignored this would be about as wide as sqrt(P[2, 2])
  conditional_sd <- stats::sd(imp$Y[5, 2, ])
  marginal_sd <- sqrt(diag(draws_mean(fit, "P"))[2])
  expect_lt(conditional_sd, 0.95 * marginal_sd)
  expect_gt(conditional_sd, 0)
})

test_that("the imputation refuses an object it cannot complete", {
  sim <- simulate_animal(n_fam = 5, size = 4, G = G3, E = E3, seed = 94)
  gibbs <- lmm_animal_gibbs(sim$Y, sim$X, sim$A, iter = 40, warmup = 20)
  expect_error(imputed_data(gibbs), "does not carry the data")
})
