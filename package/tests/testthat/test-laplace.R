# The Laplace fitters share the model code and data assembly of the samplers,
# return a summary rather than draws, and are refused by the draw-based
# functions.

skip_if_not_installed <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    testthat::skip(paste(pkg, "not installed"))
}

small_animal_data <- function(n_gen = 3, n_breed = 20) {
  set.seed(11)
  ped <- data.frame(id = seq_len(n_gen * n_breed), sire = NA_integer_,
                    dam = NA_integer_)
  for (g in 2:n_gen) {
    kids <- ((g - 1) * n_breed + 1):(g * n_breed)
    pool <- ped$id[ped$id < min(kids)]
    ped$sire[kids] <- sample(pool, n_breed, TRUE)
    ped$dam[kids]  <- sample(pool, n_breed, TRUE)
  }
  A <- makeA(ped)
  n <- nrow(A)
  G <- random_covariance(2, trait_names = c("t1", "t2"))
  E <- random_covariance(2, trait_names = c("t1", "t2"))
  X <- cbind(Intercept = 1, x = rnorm(n))
  Y <- X %*% rbind(c(0.3, -0.1), c(1, 0.4)) +
    t(chol(A)) %*% matrix(rnorm(n * 2), n) %*% chol(G) +
    matrix(rnorm(n * 2), n) %*% chol(E)
  rownames(Y) <- as.character(ped$id)
  colnames(Y) <- c("t1", "t2")
  list(Y = Y, X = X, A = A)
}

test_that("lmm_animal_laplace returns a summary and no draws to use", {
  skip_if_not_installed("cmdstanr")
  d <- small_animal_data()
  fit <- lmm_animal_laplace(d$Y, d$X, d$A, draws = 200, seed = 1)

  expect_s3_class(fit, "stanimal_laplace")
  expect_named(fit, c("summary", "draws", "elapsed", "draws_requested",
                      "seed", "call"))
  expect_true(all(c("G[1,1]", "G[2,2]", "E[1,1]", "P[1,1]", "h2[1]") %in%
                    fit$summary$quantity))
  expect_true(all(fit$summary$sd[fit$summary$quantity == "G[1,1]"] > 0))
  # the covariance is reported on the data scale, not on the standardized one
  expect_gt(fit$summary$mean[fit$summary$quantity == "P[1,1]"],
            fit$summary$mean[fit$summary$quantity == "G[1,1]"])
})

test_that("the data and model code match what the sampler would use", {
  skip_if_not_installed("cmdstanr")
  d <- small_animal_data()
  cap <- file.path(tempdir(), "stanimal_laplace_capture.rds")
  ns <- asNamespace("stanimal")

  trace("cstan", where = ns, print = FALSE,
        tracer = quote(saveRDS(list(code = model_code, data = data),
                               file.path(tempdir(),
                                         "stanimal_laplace_capture.rds"))))
  lmm_animal(d$Y, d$X, d$A, iter = 100, warmup = 50, chains = 1, cores = 1,
             refresh = 0, seed = 1)
  untrace("cstan", where = ns)
  sampler <- readRDS(cap)

  trace(".laplace_fit", where = ns, print = FALSE,
        tracer = quote(saveRDS(list(code = model_code, data = stan_data),
                               file.path(tempdir(),
                                         "stanimal_laplace_capture.rds"))))
  invisible(lmm_animal_laplace(d$Y, d$X, d$A, draws = 10, seed = 1))
  untrace(".laplace_fit", where = ns)
  laplace <- readRDS(cap)

  expect_identical(laplace$code, sampler$code)
  expect_equal(laplace$data, sampler$data)
})

test_that("univariate, factor and phenotypic Laplace fits work", {
  skip_if_not_installed("cmdstanr")
  d <- small_animal_data()

  uni <- lmm_animal_laplace(d$Y[, 1, drop = FALSE], d$X, d$A,
                            draws = 100, seed = 1)
  expect_setequal(uni$summary$quantity, c("sigma_G", "sigma_E", "h2"))

  fac <- suppressWarnings(
    lmm_animal_laplace(d$Y, d$X, d$A, genetic_structure = "factor",
                       n_factors = 1, draws = 100, seed = 1))
  expect_true(all(c("loadings[1,1]", "communality[1]", "factor_magnitude[1]") %in%
                    fac$summary$quantity))

  phen <- lm_phen_cov_laplace(d$Y, d$X, draws = 100, seed = 1)
  expect_true(all(c("P[1,1]", "corrP[2,1]", "sigma_P[1]") %in%
                    phen$summary$quantity))
})

test_that("draw-based functions refuse a Laplace object", {
  skip_if_not_installed("cmdstanr")
  d <- small_animal_data()
  fit <- lmm_animal_laplace(d$Y, d$X, d$A, draws = 100, seed = 1)

  expect_error(fit_draws(fit), "needs posterior draws")
  expect_error(draws_mean(fit, "G"), "needs posterior draws")
  expect_error(fit_loo(fit), "needs posterior draws")
  expect_error(breeding_values(fit), "does not carry the data")
})
