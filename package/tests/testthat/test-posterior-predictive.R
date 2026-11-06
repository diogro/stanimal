# The tail probability is the share of replicates at or above the observation,
# which is a property of the definition and can be checked without a sampler: an
# observation sitting at the centre of the replicates must land near 0.5, and one
# outside them must land at the ends.
test_that("the predictive tail probability counts replicates above the data", {
  expect_equal(stanimal:::.predictive_tail(c(1, 2, 3, 4, 5), 3), 0.6)
  expect_equal(stanimal:::.predictive_tail(c(1, 2, 3), 10), 0)
  expect_equal(stanimal:::.predictive_tail(c(1, 2, 3), 0), 1)

  set.seed(1)
  replicated <- stats::rnorm(1000)
  expect_equal(stanimal:::.predictive_tail(replicated, stats::median(replicated)),
               0.5, tolerance = 0.05)
  # a feature the model cannot produce is at the ends, not in the middle
  expect_lt(stanimal:::.predictive_tail(replicated,
                                        stats::quantile(replicated, 0.999)), 0.01)
})

predictive_mock_fit <- function(A, G, E, beta, ndraws) {
  k <- nrow(G)
  n <- nrow(A)
  names_matrix <- function(name) {
    sprintf("%s[%d,%d]", name, rep(seq_len(k), k), rep(seq_len(k), each = k))
  }
  values <- c(as.vector(G), as.vector(E), beta)
  vars <- c(names_matrix("G"), names_matrix("E"),
            sprintf("beta[%d,1]", seq_len(k)))
  m <- matrix(rep(values, each = ndraws), ndraws,
              dimnames = list(NULL, vars))
  Y <- matrix(0, n, k, dimnames = list(paste0("id", seq_len(n)),
                                       paste0("trait", seq_len(k))))
  structure(list(draws = function(...) m),
            class = "CmdStanMCMC", y_names = colnames(Y), x_names = "Int",
            stanimal_data = list(family = "animal", Y = Y,
                                 X = matrix(1, n, 1), A = A))
}

test_that("animal replicates use a left factor of the relationship matrix", {
  draws <- list(matrix(c(0.2, -0.3), 2), matrix(c(1, -2), 2))
  calls <- 0L
  local_mocked_bindings(.rmvn_chol = function(n, k, Sigma) {
    calls <<- calls + 1L
    draws[[calls]]
  }, .package = "stanimal")
  for (A in list(matrix(c(1.5, 0.6, 0.6, 0.8), 2), diag(2))) {
    calls <- 0L
    fit <- predictive_mock_fit(A, matrix(2), matrix(0.5), 3, 1)
    result <- posterior_predictive(fit, ndraws = 1, seed = 1)
    expected <- 3 + draws[[1]] + t(chol(A + stanimal:::.jitter(A))) %*% draws[[2]]
    expect_equal(unname(result$yrep[, , 1]), as.numeric(expected))
    expect_equal(dim(result$yrep), c(2L, 1L, 1L))
    expect_equal(dimnames(result$yrep)[1:2], dimnames(attr(fit, "stanimal_data")$Y))
  }
})

test_that("univariate SD draws and phenotypic fits retain their simulation paths", {
  A <- matrix(c(1.5, 0.6, 0.6, 0.8), 2)
  Y <- matrix(0, 2, 1)
  data <- list(family = "animal", Y = Y, X = matrix(1, 2, 1), A = A)
  draws <- matrix(c(2, 0.5, 3), 1,
                  dimnames = list(NULL, c("sigma_G", "sigma_E", "beta[1]")))
  fit <- structure(list(draws = function(...) draws), stanimal_data = data)
  set.seed(87)
  residual <- stanimal:::.rmvn_chol(2, 1, matrix(0.25))
  genetic <- stanimal:::.rmvn_chol(2, 1, matrix(4))
  expected <- 3 + residual + t(chol(A + stanimal:::.jitter(A))) %*% genetic
  expect_equal(as.numeric(posterior_predictive(fit, 1, seed = 87)$yrep),
               as.numeric(expected))

  attr(fit, "stanimal_data")$family <- "phenotypic"
  draws <- matrix(c(2, 3), 1, dimnames = list(NULL, c("P[1,1]", "beta[1]")))
  fit$draws <- function(...) draws
  set.seed(87)
  expected <- 3 + stanimal:::.rmvn_chol(2, 1, matrix(2))
  expect_equal(as.numeric(posterior_predictive(fit, 1, seed = 87)$yrep),
               as.numeric(expected))
})

test_that("animal replicate covariance matches A tensor G plus residual E", {
  A <- matrix(c(1.5, 0.6, 0.2, 0.6, 0.8, 0.1, 0.2, 0.1, 1.2), 3)
  G <- matrix(c(1.2, 0.35, 0.35, 0.7), 2)
  E <- matrix(c(0.4, -0.12, -0.12, 0.6), 2)
  fit <- predictive_mock_fit(A, G, E, c(2, -1), 20000)
  pp <- posterior_predictive(fit, ndraws = 20000, seed = 271)
  empirical <- stats::cov(t(matrix(pp$yrep, nrow = 6)))
  expected <- kronecker(G, A) + kronecker(E, diag(3))
  # The fixed-seed maximum deviation is 0.045; allow twice that Monte Carlo error.
  expect_lt(max(abs(empirical - expected)), 0.09)
  expect_equal(dim(pp$yrep), c(3L, 2L, 20000L))
  expect_equal(pp$yrep, posterior_predictive(fit, ndraws = 20000, seed = 271)$yrep)
})

test_that("predictive statistics use the observed mask for each trait and pair", {
  Y <- cbind(a = c(1, 2, 3, NA, NA), b = c(NA, 3, 5, 6, NA))
  replicate <- cbind(a = c(1, 2, 3, 20, 20), b = c(99, 3, 5, 6, 99))
  x <- structure(list(y = Y, yrep = array(rep(replicate, 2), c(5, 2, 2)),
                      traits = c("a", "b"), ndraws = 2L),
                 class = "stanimal_yrep")
  s <- stanimal:::.yrep_statistics(x)
  for (name in unique(s$statistic)) {
    item <- s[s$statistic == name, ]
    expect_equal(item$replicated, rep(item$observed[1], 2))
  }
  expect_equal(s$observed[s$statistic == "mean: a"][1], 2)
  expect_output(print(x), "posterior predictive")
  expect_s3_class(plot(x), "ggplot")
})

test_that("unavailable pairwise correlations are reported, not silently omitted", {
  Y <- cbind(a = c(1, 2, NA, NA), b = c(NA, NA, 3, 4))
  x <- structure(list(y = Y, yrep = array(seq_len(16), c(4, 2, 2)),
                      traits = c("a", "b"), ndraws = 2L),
                 class = "stanimal_yrep")
  expect_warning(s <- stanimal:::.yrep_statistics(x), "correlation.*a.*b")
  expect_false(any(grepl("^cor:", s$statistic)))
  x$y <- cbind(a = c(1, 1, 2, NA), b = c(3, 4, NA, 5))
  expect_warning(s <- stanimal:::.yrep_statistics(x), "correlation.*a.*b")
  expect_false(any(grepl("^cor:", s$statistic)))
  x$y <- cbind(a = 1:4, b = 4:1)
  x$yrep[, 1, 1] <- 1
  expect_error(stanimal:::.yrep_statistics(x), "replicated correlation.*a-b")
})

test_that("a well-specified fit does not flag the features it reproduces", {
  skip_if_no_cmdstan()
  sim <- simulate_animal(n_fam = 12, size = 7, G = G3, E = E3, seed = 81)
  fit <- suppressWarnings(lmm_animal(sim$Y, sim$X, sim$A, iter = 600,
                                     chains = 2, cores = 2, seed = 1,
                                     refresh = 0, show_messages = FALSE))
  pp <- posterior_predictive(fit, ndraws = 100, seed = 2)
  out <- capture.output(print(pp))
  # the table's p column: parse it rather than re-deriving it, so the test
  # covers what the user reads
  rows <- grep("^ +(mean|sd|cor): ", out, value = TRUE)
  expect_length(rows, 3L + 3L + 3L)
  p <- as.numeric(sub("^.* ([0-9.]+)$", "\\1", rows))
  # This particular simulated example should not flag most features.
  expect_equal(sum(p == 0 | p == 1), 0L)
  expect_gte(sum(p > 0.1 & p < 0.9), 5L)
})
