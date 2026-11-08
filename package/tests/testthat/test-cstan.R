test_that("cstan fits a Stan program and returns a CmdStanMCMC", {
  skip_if_no_cmdstan()
  code <- "
  data { int<lower=0> N; vector[N] y; }
  parameters { real mu; real<lower=0> sigma; }
  model { y ~ normal(mu, sigma); mu ~ normal(0, 10); sigma ~ normal(0, 5); }"
  set.seed(1)
  y <- stats::rnorm(50, 3, 2)
  fit <- stanimal:::cstan(code, list(N = 50, y = y), chains = 2, cores = 2,
                            iter = 800, refresh = 0, show_messages = FALSE,
                            seed = 1)
  expect_s3_class(fit, "CmdStanMCMC")
  expect_equal(posterior::ndraws(fit$draws()), 800)
  expect_lt(abs(fit$summary("mu", "mean")$mean - 3), 0.8)
})

test_that("cstan passes adapt_delta and max_treedepth from control", {
  skip_if_no_cmdstan()
  code <- "parameters { real x; } model { x ~ normal(0, 1); }"
  fit <- stanimal:::cstan(code, list(), chains = 1, iter = 200,
                            control = list(adapt_delta = 0.8, max_treedepth = 5),
                            refresh = 0, show_messages = FALSE)
  expect_equal(fit$metadata()$adapt_delta, 0.8)
  expect_equal(fit$metadata()$max_treedepth, 5)
})

test_that("chains start from init = 1 unless init is passed", {
  skip_if_no_cmdstan()
  # The engines' default U(-2, 2) stalls correlation factors over many traits.
  code <- "parameters { real x; } model { x ~ normal(0, 1); }"
  fit <- function(...) suppressWarnings(stanimal:::cstan(
    code, list(), chains = 1, iter = 100, refresh = 0, show_messages = FALSE, ...))
  expect_equal(fit()$metadata()$init, 1)
  expect_equal(fit(init = 0.5)$metadata()$init, 0.5)
})

test_that("cstan forwards the other sampler controls", {
  skip_if_no_cmdstan()
  code <- "parameters { real x; } model { x ~ normal(0, 1); }"
  fit <- stanimal:::cstan(code, chains = 1, iter = 200, refresh = 0,
                            show_messages = FALSE,
                            control = list(adapt_engaged = FALSE, step_size = 0.3))
  stepsize <- as.matrix(fit$sampler_diagnostics(format = "draws_matrix"))[, "stepsize__"]
  expect_equal(as.numeric(unique(stepsize)), 0.3)
  expect_error(stanimal:::cstan(code, iter = 100, control = list(gamma = 0.5)),
               "unsupported sampler controls: gamma")
})

test_that("control is validated once, for both engines", {
  expect_equal(stanimal:::.check_control(NULL), list())
  expect_error(stanimal:::.check_control(list(gamma = 1)), "gamma")
  expect_error(stanimal:::.check_control(list(1, 2)), "named list")
  expect_named(stanimal:::.rstan_control(list(step_size = 0.3, adapt_delta = 0.9)),
               c("stepsize", "adapt_delta"))
  expect_named(stanimal:::.check_control(list(stepsize = 0.3)), "step_size")
  expect_error(stanimal:::.check_control(list(step_size = 0.3, stepsize = 0.4)),
               "more than once")
})

test_that("fits from the same random state keep their own output files", {
  skip_if_no_cmdstan()
  # cmdstanr names output from the minute, the pid and one runif(): without a
  # unique basename these two fits share a CSV, fit 1 reads fit 2's draws, and
  # collecting fit 1 deletes fit 2's output.
  A <- make_A(10, 5)
  set.seed(10)
  Y1 <- t(chol(A)) %*% matrix(rnorm(100), 50) + matrix(rnorm(100), 50)
  Y2 <- Y1 + 5
  fit <- function(Y) {
    set.seed(1)
    lmm_animal(Y, NULL, A, chains = 1, iter = 200, seed = 3, init = 0.5,
               refresh = 0, show_messages = FALSE)
  }
  f1 <- suppressWarnings(fit(Y1))
  f2 <- suppressWarnings(fit(Y2))

  expect_false(any(f1$output_files() %in% f2$output_files()))
  beta <- function(f) colMeans(posterior::as_draws_matrix(f$draws("beta")))
  expect_equal(unname(beta(f1)), unname(colMeans(Y1)), tolerance = 0.5)
  expect_equal(unname(beta(f2)), unname(colMeans(Y2)), tolerance = 0.5)

  files2 <- f2$output_files()
  rm(f1); invisible(gc())
  expect_true(all(file.exists(files2)))
})

test_that("a caller's output location is kept", {
  args <- stanimal:::.unique_output(list(output_dir = "somewhere"), "m.stan")
  expect_null(args$output_basename)
  args <- stanimal:::.unique_output(list(output_basename = "mine"), "m.stan")
  expect_identical(args$output_basename, "mine")
  set.seed(1); a <- stanimal:::.unique_output(list(), "m.stan")$output_basename
  set.seed(1); b <- stanimal:::.unique_output(list(), "m.stan")$output_basename
  expect_false(identical(a, b))
  expect_match(a, "^m-")
})
