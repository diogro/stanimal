# The Laplace fitters share the model code and data assembly of the samplers,
# return a summary rather than draws, and are refused by the draw-based
# functions.

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
  skip_if_no_cmdstan()
  d <- small_animal_data()
  fit <- lmm_animal_laplace(d$Y, d$X, d$A, draws = 200, seed = 1)

  expect_s3_class(fit, "stanimal_laplace")
  expect_named(fit, c("summary", "draws", "y_names", "x_names", "elapsed",
                      "draws_requested", "seed", "call"))
  expect_named(fit$summary, c("quantity", "mean", "sd", "mad"))
  expect_equal(fit$seed, 1)
  expect_equal(fit$draws_requested, 200)
  expect_equal(dim(fit$draws)[1], 200)
  # quantities are labelled with the trait names of Y (t1, t2)
  expect_true(all(c("G[t1,t1]", "G[t2,t2]", "E[t1,t1]", "P[t1,t1]",
                    "h2[t1]") %in% fit$summary$quantity))
  expect_true(all(fit$summary$sd[fit$summary$quantity == "G[t1,t1]"] > 0))
  # the covariance is reported on the data scale, not on the standardized one
  expect_gt(fit$summary$mean[fit$summary$quantity == "P[t1,t1]"],
            fit$summary$mean[fit$summary$quantity == "G[t1,t1]"])
})

test_that("the data and model code match what the sampler would use", {
  skip_if_no_cmdstan()
  d <- small_animal_data()
  seen <- new.env(parent = emptyenv())

  # Both engines are stubbed out, so nothing is compiled or sampled. What this
  # test checks is the pair of arguments each entry point hands to the Stan
  # engine, and running a real chain would make the result depend on whether
  # CmdStan's initialization happens to survive its LKJ Cholesky check. With one
  # or two chains that check failing takes the whole fit down, and the report
  # would be a missing output CSV rather than anything about the arguments.
  local_mocked_bindings(
    cstan = function(model_code, data, ...) {
      seen$sampler <- list(code = model_code, data = data)
      stop("stanimal_capture")
    },
    .package = "stanimal"
  )
  expect_error(lmm_animal(d$Y, d$X, d$A), "stanimal_capture")

  local_mocked_bindings(
    .laplace_fit = function(model_code, stan_data, ...) {
      seen$laplace <- list(code = model_code, data = stan_data)
      stop("stanimal_capture")
    },
    .package = "stanimal"
  )
  expect_error(lmm_animal_laplace(d$Y, d$X, d$A), "stanimal_capture")

  expect_identical(seen$laplace$code, seen$sampler$code)
  expect_equal(seen$laplace$data, seen$sampler$data)
})

test_that("univariate, factor and phenotypic Laplace fits work", {
  skip_if_no_cmdstan()
  d <- small_animal_data()

  uni <- lmm_animal_laplace(d$Y[, 1, drop = FALSE], d$X, d$A,
                            draws = 100, seed = 1)
  expect_setequal(uni$summary$quantity, c("sigma_G", "sigma_E", "h2"))

  fac <- lmm_animal_laplace(d$Y, d$X, d$A, genetic_structure = "factor",
                            n_factors = 1, draws = 100, seed = 1)
  expect_true(all(c("loadings[t1,1]", "communality[t1]", "factor_magnitude[1]") %in%
                    fac$summary$quantity))

  phen <- lm_phen_cov_laplace(d$Y, d$X, draws = 100, seed = 1)
  expect_true(all(c("P[t1,t1]", "corrP[t2,t1]", "sigma_P[t1]") %in%
                    phen$summary$quantity))
  phen_again <- lm_phen_cov_laplace(d$Y, d$X, draws = 100, seed = 1)
  expect_equal(phen$summary, phen_again$summary)
})

test_that("draw-based functions refuse a Laplace object", {
  skip_if_no_cmdstan()
  d <- small_animal_data()
  fit <- lmm_animal_laplace(d$Y, d$X, d$A, draws = 100, seed = 1)

  expect_error(fit_draws(fit), "needs posterior draws")
  expect_error(draws_mean(fit, "G"), "needs posterior draws")
  expect_error(fit_loo(fit), "needs posterior draws")
  expect_error(breeding_values(fit), "does not carry the data")
  expect_error(posterior_predictive(fit), "does not carry the data")
  expect_error(imputed_data(fit), "does not carry the data")
  expect_error(plot_prior_posterior(fit), "needs posterior draws")
  expect_error(plot_forest(fit), "needs posterior draws")
})

test_that("unused arguments to the Laplace fitters warn", {
  expect_silent(stanimal:::.warn_dots())
  expect_warning(stanimal:::.warn_dots(iter = 10), "iter")
  expect_warning(stanimal:::.warn_dots(10, iter = 10), "<unnamed>, iter")
})

test_that("Laplace failures and warnings are surfaced without guessing", {
  calls <- 0L
  model <- list(laplace = function(...) {
    calls <<- calls + 1L
    warning("optimizer warning")
    stop("unrelated optimizer error")
  })
  expect_warning(expect_error(stanimal:::.laplace_run(model, list(), "h2", 10, 1),
                              "unrelated optimizer error"), "optimizer warning")
  expect_equal(calls, 1L)

  model$laplace <- function(...) list(
    return_codes = function() 1L,
    output = function(...) cat("diagnostic: unit vector is invalid"))
  expect_error(stanimal:::.laplace_run(model, list(), "h2", 10, 1),
               "diagnostic: unit vector is invalid")

  model$laplace <- function(...) {
    args <- list(...)
    expect_equal(args$seed, 9)
    expect_equal(args$sig_figs, 18)
    warning("successful optimizer warning")
    list(return_codes = function() 0L,
         summary = function(...) data.frame(variable = "h2", mean = 0.4,
                                             sd = 0.1, mad = 0.1),
         draws = function(...) matrix(0.4))
  }
  expect_warning(result <- stanimal:::.laplace_run(model, list(), "h2", 10, 9),
                 "successful optimizer warning")
  expect_equal(result$seed, 9)
  expect_equal(result$summary$mean, 0.4)
})

test_that("both Laplace fitters match sampler model code and data", {
  d <- small_animal_data()
  d$Y[1, 1] <- NA_real_
  seen <- new.env(parent = emptyenv())
  local_mocked_bindings(
    .resolve_stan_backend = function(cmd_stan) TRUE,
    cstan = function(model_code, data, ...) {
      seen$sampler <- list(code = model_code, data = data)
      stop("capture sampler")
    },
    .laplace_fit = function(model_code, stan_data, ...) {
      seen$laplace <- list(code = model_code, data = stan_data)
      stop("capture Laplace")
    },
    .package = "stanimal"
  )
  cases <- list(
    list(Y = d$Y[, 1, drop = FALSE], args = list()),
    list(Y = d$Y[, 1, drop = FALSE], args = list(parametrization = "h2",
                                                  beta_prior = c(2, 5))),
    list(Y = d$Y, args = list()),
    list(Y = d$Y, args = list(parametrization = "h2", beta_prior = c(2, 5),
                               residual_structure = "diagonal")),
    list(Y = d$Y, args = list(genetic_structure = "factor", n_factors = 1)),
    list(Y = d$Y, args = list(genetic_structure = "factor", n_factors = 1,
                               residual_structure = "diagonal",
                               hs_global_scale = 0.7, hs_slab_scale = 1.5,
                               column_decay = 0.6, beta_prior = c(2, 3))),
    list(Y = d$Y, args = list(genetic_structure = "factor", n_factors = 1,
                               residual_structure = "shared_factor"))
  )
  for (case in cases) {
    args <- c(list(Y = case$Y, X = d$X, A = d$A), case$args)
    expect_error(do.call(lmm_animal, args), "capture sampler")
    expect_error(do.call(lmm_animal_laplace, args), "capture Laplace")
    expect_identical(seen$laplace$code, seen$sampler$code)
    expect_equal(seen$laplace$data, seen$sampler$data)
  }
  for (structure in c("unstructured", "factor")) {
    args <- list(Y = d$Y, X = d$X, covariance_structure = structure)
    if (structure == "factor") args$n_factors <- 1
    expect_error(do.call(lm_phen_cov, args), "capture sampler")
    expect_error(do.call(lm_phen_cov_laplace, args), "capture Laplace")
    expect_identical(seen$laplace$code, seen$sampler$code)
    expect_equal(seen$laplace$data, seen$sampler$data)
  }
})

test_that("Laplace fitters warn about unused model arguments", {
  d <- small_animal_data()
  local_mocked_bindings(.laplace_fit = function(...) stop("capture"),
                        .package = "stanimal")
  expect_warning(expect_error(lmm_animal_laplace(d$Y, d$X, d$A,
                                                   beta_prior = c(2, 2)),
                              "capture"), "`beta_prior`")
  expect_no_warning(expect_error(lmm_animal_laplace(
    d$Y, d$X, d$A, residual_structure = "diagonal", lkj_prior = 3),
    "capture"))
  expect_warning(expect_error(lm_phen_cov_laplace(
    d$Y, d$X, covariance_structure = "factor", n_factors = 1, lkj_prior = 3),
    "capture"), "`lkj_prior`")
  expect_warning(expect_error(lmm_animal_laplace(
    d$Y, d$X, d$A, genetic_structure = "factor", n_factors = 1,
    residual_structure = "shared_factor", lkj_prior = 3), "capture"),
    "`lkj_prior`")
  expect_warning(expect_error(lmm_animal_laplace(
    d$Y, d$X, d$A, genetic_structure = "factor", n_factors = 1,
    parametrization = "h2"), "capture"), "`parametrization`")
})

test_that("LKJ use reflects both genetic and residual structures in both fitters", {
  d <- small_animal_data()
  local_mocked_bindings(
    .resolve_stan_backend = function(cmd_stan) TRUE,
    cstan = function(...) stop("capture sampler"),
    .laplace_fit = function(...) stop("capture Laplace"),
    .package = "stanimal"
  )
  cases <- list(
    list(Y = d$Y[, 1, drop = FALSE], args = list(), used = FALSE),
    list(Y = d$Y, args = list(residual_structure = "unstructured"), used = TRUE),
    list(Y = d$Y, args = list(residual_structure = "diagonal"), used = TRUE),
    list(Y = d$Y, args = list(genetic_structure = "factor", n_factors = 1,
                              residual_structure = "unstructured"), used = TRUE),
    list(Y = d$Y, args = list(genetic_structure = "factor", n_factors = 1,
                              residual_structure = "diagonal"), used = FALSE),
    list(Y = d$Y, args = list(genetic_structure = "factor", n_factors = 1,
                              residual_structure = "shared_factor"), used = FALSE)
  )
  for (case in cases) for (fn in list(lmm_animal, lmm_animal_laplace)) {
    call <- function() do.call(fn, c(list(Y = case$Y, X = d$X, A = d$A,
                                          lkj_prior = 3), case$args))
    if (case$used) expect_no_warning(expect_error(call(), "capture"))
    else expect_warning(expect_error(call(), "capture"), "`lkj_prior`")
  }
})

test_that("Laplace validates configurations and run arguments", {
  d <- small_animal_data()
  expect_error(lmm_animal_laplace(d$Y, d$X, d$A,
                                   residual_structure = "shared_factor"),
               "requires")
  expect_error(lmm_animal_laplace(d$Y, d$X, d$A, n_factors = 1),
               "only used")
  expect_error(lm_phen_cov_laplace(d$Y, d$X,
                                   covariance_structure = "factor",
                                   n_factors = 2), "smaller")
  expect_error(lm_phen_cov_laplace(d$Y, d$X, lkj_prior = 0), "lkj_prior")
  expect_error(lm_phen_cov_laplace(d$Y, d$X, draws = 0), "draws")
  expect_error(lm_phen_cov_laplace(d$Y, d$X, seed = -1), "seed")
})

test_that("additional Laplace model structures run with missing responses", {
  skip_if_no_cmdstan()
  d <- small_animal_data()
  d$Y[1, 1] <- NA_real_
  uni <- lmm_animal_laplace(d$Y[, 1, drop = FALSE], d$X, d$A,
                            parametrization = "h2", draws = 100, seed = 1)
  expect_true("h2" %in% uni$summary$quantity)
  diagonal <- lmm_animal_laplace(d$Y, d$X, d$A,
                                 residual_structure = "diagonal",
                                 draws = 100, seed = 1)
  expect_true("G[t1,t1]" %in% diagonal$summary$quantity)
  shared <- lmm_animal_laplace(d$Y, d$X, d$A,
                               genetic_structure = "factor", n_factors = 1,
                               residual_structure = "shared_factor",
                               draws = 100, seed = 1)
  expect_true("factor_heritability[1]" %in% shared$summary$quantity)
  phen <- lm_phen_cov_laplace(d$Y, d$X, covariance_structure = "factor",
                              n_factors = 1, draws = 100, seed = 1)
  expect_true("P[t1,t1]" %in% phen$summary$quantity)
})

test_that("Laplace and sampling locate the same region on a regular example", {
  skip_if_no_cmdstan()
  set.seed(2026)
  n <- 200
  X <- cbind(intercept = 1, predictor = rnorm(n))
  P <- matrix(c(1, 0.4, 0.4, 1.5), 2)
  Y <- X %*% rbind(c(0.5, -0.5), c(0.3, -0.2)) +
    matrix(rnorm(n * 2), n) %*% chol(P)
  approx <- lm_phen_cov_laplace(Y, X, draws = 2000, seed = 41)
  sampled <- lm_phen_cov(Y, X, iter = 1000, warmup = 500, chains = 2,
                         cores = 2, seed = 41, refresh = 0)
  quantities <- c("P[1,1]", "P[2,2]", "corrP[2,1]")
  a <- approx$summary[match(quantities, approx$summary$quantity), ]
  s <- sampled$summary(quantities)
  s <- s[match(quantities, s$variable), ]
  expect_true(all(s$rhat < 1.05))
  expect_equal(sum(sampled$diagnostic_summary()$num_divergent), 0)
  # On this fixed example the measured relative diagonal differences are
  # below 0.03 and the correlation difference is below 0.004.
  expect_lt(max(abs(a$mean[1:2] / s$mean[1:2] - 1)), 0.06)
  expect_lt(abs(a$mean[3] - s$mean[3]), 0.01)
})

test_that("Laplace quantities are labelled with the trait and predictor names", {
  label <- stanimal:::.label_quantities
  y <- c("tail", "wing"); x <- c("Int", "size")
  expect_equal(label(c("G[2,1]", "corrP[1,2]", "h2[2]", "sigma_G", "beta[2,1]",
                       "loadings[1,3]", "lp__", "factor_magnitude[1]"), y, x),
               c("G[wing,tail]", "corrP[tail,wing]", "h2[wing]", "sigma_G",
                 "beta[wing,Int]", "loadings[tail,3]", "lp__",
                 "factor_magnitude[1]"))
  # no names: unchanged; an index beyond the names: unchanged
  expect_equal(label(c("G[2,1]", "h2[1]"), NULL), c("G[2,1]", "h2[1]"))
  expect_equal(label("G[3,1]", y), "G[3,1]")

  skip_if_no_cmdstan()
  d <- simulate_animal(n_fam = 8, size = 5, G = diag(2), E = diag(2))
  Y <- d$Y
  colnames(Y) <- c("tail", "wing")
  fit <- lmm_animal_laplace(Y, d$X, d$A, draws = 100, seed = 1)
  expect_true(all(c("G[tail,tail]", "G[wing,tail]", "h2[wing]") %in%
                    fit$summary$quantity))
  expect_true("G[wing,tail]" %in% posterior::variables(fit$draws))
  expect_identical(fit$y_names, c("tail", "wing"))
})
