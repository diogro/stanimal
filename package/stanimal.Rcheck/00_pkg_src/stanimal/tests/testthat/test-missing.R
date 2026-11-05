# Incomplete trait records.
#
# The phenotypic models integrate a pattern's missing values out analytically;
# the animal models cannot (a gap is informed by every relative), so they fill
# the gaps with parameters. Both must reproduce the likelihood of the observed
# data exactly, which is what the density checks below verify.

# Everything the missing-response helper promises, without touching Stan.
test_that("missing values are grouped into patterns", {
  Y <- matrix(c(1, 2, 3,
                4, NA, 6,
                7, 8, NA,
                10, 11, 12,
                13, NA, 15), 5, 3, byrow = TRUE)
  m <- stanimal:::.missing_response(Y)

  expect_equal(m$n_missing, 3L)
  expect_equal(m$n_pattern, 3L)          # {1,2,3}, {1,3}, {1,2}
  expect_equal(m$n_rows, c(2L, 2L, 1L))  # rows 1 and 4 are complete
  expect_equal(m$pattern_obs, matrix(c(1L, 1L, 1L,
                                       1L, 0L, 1L,
                                       1L, 1L, 0L), 3, 3, byrow = TRUE))
  # rows are reordered by pattern, so every pattern is one contiguous block
  expect_equal(m$row_order, c(1L, 4L, 2L, 5L, 3L))
  expect_equal(m$first_row, c(1L, 3L, 5L))
  # no NA reaches Stan, and the observed entries are untouched
  expect_false(anyNA(m$Y))
  expect_equal(m$Y[!is.na(Y)], Y[!is.na(Y)])
  expect_equal(m$Y[is.na(Y)], c(0, 0, 0))
  # the gaps are numbered in the order the parameters are declared, which is
  # column-major, the same order as the Stan vector
  expect_equal(m$mis_row, c(2L, 5L, 3L))
  expect_equal(m$mis_col, c(2L, 2L, 3L))
})

test_that("a response with no gaps is a single complete pattern", {
  Y <- matrix(stats::rnorm(12), 4, 3)
  m <- stanimal:::.missing_response(Y)
  expect_equal(m$n_missing, 0L)
  expect_equal(m$n_pattern, 1L)
  expect_equal(m$n_rows, 4L)
  expect_equal(m$row_order, 1:4)
  expect_equal(m$pattern_obs, matrix(1L, 1, 3))
})

test_that("unusable responses are refused with the reason", {
  Y <- matrix(stats::rnorm(12), 4, 3)
  Y[, 2] <- NA
  expect_error(stanimal:::.missing_response(Y), "at least one observed value")
  Y <- matrix(stats::rnorm(12), 4, 3)
  Y[2, ] <- NA
  expect_error(stanimal:::.missing_response(Y), "at least one observed trait")
  # the animal models can take such a record, because relatives inform it
  expect_no_error(stanimal:::.missing_response(Y, allow_empty_rows = TRUE))
})

test_that("animal_design keeps records that miss only the response", {
  set.seed(4)
  dat <- data.frame(id = as.character(1:6), y1 = stats::rnorm(6),
                    y2 = stats::rnorm(6), sex = factor(c("f", "m")))
  dat$y2[2] <- NA
  dat$sex[4] <- NA
  rownames(dat) <- dat$id
  design <- animal_design(cbind(y1, y2) ~ sex, dat)
  # the row with an unknown predictor goes, the one with a missing trait stays
  expect_equal(rownames(design$Y), c("1", "2", "3", "5", "6"))
  expect_true(is.na(design$Y["2", "y2"]))
  expect_equal(nrow(design$X), nrow(design$Y))
})

test_that("phenotypic fits use incomplete records instead of dropping them", {
  skip_if_no_cmdstan()
  set.seed(5)
  n <- 150
  k <- 3
  sds <- c(1, 2, 3)
  Rmat <- matrix(0.5, k, k)
  diag(Rmat) <- 1
  P <- diag(sds) %*% Rmat %*% diag(sds)
  X <- cbind(Int = 1, x = stats::rnorm(n))
  B <- cbind(c(1, 2, 3), c(-1, 0, 1))
  Y <- X %*% t(B) + matrix(stats::rnorm(n * k), n, k) %*% chol(P)
  colnames(Y) <- c("t1", "t2", "t3")

  fit_complete <- lm_phen_cov(Y, X, iter = 800, chains = 2, cores = 2, seed = 1,
                              refresh = 0, show_messages = FALSE)
  set.seed(6)
  Y_miss <- Y
  for (i in 1:40) {                       # whole rows of one trait at a time
    j <- sample(k, 1)
    if (sum(!is.na(Y_miss[i, ])) > 1) Y_miss[i, j] <- NA
  }
  expect_gt(sum(is.na(Y_miss)), 20)

  fit_missing <- lm_phen_cov(Y_miss, X, iter = 800, chains = 2, cores = 2,
                             seed = 1, refresh = 0, show_messages = FALSE)
  P_complete <- unname(draws_mean(fit_complete, "P"))
  P_missing <- unname(draws_mean(fit_missing, "P"))
  expect_lt(max(abs(P_missing - P_complete)) / max(abs(P_complete)), 0.15)
  expect_equal(dim(fit_missing$draws("log_lik", format = "draws_matrix")),
               c(800L, n))
})

test_that("the phenotypic density is the observed-data marginal", {
  skip_if_no_cmdstan()
  set.seed(8)
  n <- 60
  k <- 3
  Rmat <- matrix(0.4, k, k)
  diag(Rmat) <- 1
  P <- diag(c(1, 2, 3)) %*% Rmat %*% diag(c(1, 2, 3))
  Y <- cbind(c(1, 2, 3), c(1.5, 2, 2.5), c(0.5, 2, 3.5))[rep(1:3, length.out = n), ] +
    matrix(stats::rnorm(n * k), n, k) %*% chol(P)
  Y[c(3, 9, 14), 2] <- NA
  Y[c(1, 7), 3] <- NA
  Y[5, c(1, 3)] <- NA

  # intercept only, so the model's design is the standardized constant column
  fit <- lm_phen_cov(Y, NULL, iter = 400, chains = 1, seed = 2, refresh = 0,
                     show_messages = FALSE)
  d <- 5L
  get <- function(v) as.numeric(fit$draws(v, format = "draws_matrix")[d, ])
  beta_scaled <- get("beta_scaled")
  corrP <- matrix(0, k, k)
  for (i in 1:k) for (j in 1:k) corrP[i, j] <- get(sprintf("corrP[%d,%d]", i, j))
  y_sd <- apply(Y, 2L, stats::sd, na.rm = TRUE)
  s_P <- get("sigma_P") / y_sd
  P_std <- diag(s_P) %*% corrP %*% diag(s_P)

  # the model block works on the standardized response, so the comparison is
  # there; the raw-scale density differs from it by a constant Jacobian
  Y_std <- Y
  for (j in 1:k) {
    obs <- !is.na(Y[, j])
    Y_std[, j] <- (Y[, j] - mean(Y[obs, j])) / stats::sd(Y[obs, j])
  }
  resid <- Y_std - rep(beta_scaled, each = n)

  key <- apply(!is.na(Y), 1L, function(o) paste(as.integer(o), collapse = ""))
  expected <- 0
  for (pattern in unique(key)) {
    rows <- which(key == pattern)
    obs <- which(!is.na(Y[rows[1], ]))
    L <- chol(P_std[obs, obs, drop = FALSE])
    z <- forwardsolve(t(L), t(resid[rows, obs, drop = FALSE]))
    expected <- expected + sum(stats::dnorm(z, log = TRUE)) -
      length(rows) * sum(log(diag(L)))
  }
  expect_equal(sum(get("log_lik")), expected, tolerance = 1e-6)
})

test_that("animal fits impute gaps from relatives", {
  skip_if_no_cmdstan()
  sim <- simulate_animal(n_fam = 12, size = 6, G = G3, E = E3, seed = 12)
  A <- sim$A
  X <- sim$X

  fit_complete <- lmm_animal(sim$Y, X, A, iter = 800, chains = 2, cores = 2,
                             seed = 1, refresh = 0, show_messages = FALSE)
  P_complete <- unname(draws_mean(fit_complete, "P"))

  set.seed(13)
  Y_miss <- sim$Y
  for (i in 1:30) {
    j <- sample(ncol(Y_miss), 1)
    if (sum(!is.na(Y_miss[i, ])) > 1) Y_miss[i, j] <- NA
  }
  fit_missing <- lmm_animal(Y_miss, X, A, iter = 800, chains = 2, cores = 2,
                            seed = 1, refresh = 0, show_messages = FALSE)
  P_missing <- unname(draws_mean(fit_missing, "P"))
  expect_lt(max(abs(P_missing - P_complete)) / max(abs(P_complete)), 0.3)

  # the gap parameters are reported, so the imputed values can be inspected
  expect_true("Y_mis" %in% fit_missing$metadata()$stan_variables)
  expect_equal(dim(fit_missing$draws("Y_mis", format = "draws_matrix")),
               c(800L, sum(is.na(Y_miss))))
})

test_that("the animal density is the joint density of the completed data", {
  skip_if_no_cmdstan()
  sim <- simulate_animal(n_fam = 8, size = 5, G = G3, E = E3, seed = 14)
  A <- sim$A
  Y <- sim$Y
  Y[c(2, 11), 1] <- NA
  Y[c(7), 3] <- NA
  Y[20, ] <- NA                          # nothing observed: imputed from relatives

  fit <- lmm_animal(Y, NULL, A, iter = 300, chains = 1, seed = 3, refresh = 0,
                    show_messages = FALSE)
  d <- 4L
  get <- function(v) as.numeric(fit$draws(v, format = "draws_matrix")[d, ])
  k <- ncol(Y)
  n <- nrow(Y)
  L_G <- diag(k); L_E <- diag(k)
  for (i in 1:k) for (j in 1:i) {
    L_G[i, j] <- get(sprintf("L_Omega_G[%d,%d]", i, j))
    L_E[i, j] <- get(sprintf("L_Omega_E[%d,%d]", i, j))
  }
  G_std <- diag(get("s_G")) %*% tcrossprod(L_G) %*% diag(get("s_G"))
  E_std <- diag(get("s_E")) %*% tcrossprod(L_E) %*% diag(get("s_E"))

  miss <- is.na(Y)
  Y_std <- Y
  for (j in 1:k) {
    obs <- !miss[, j]
    Y_std[, j] <- (Y[, j] - mean(Y[obs, j])) / stats::sd(Y[obs, j])
  }
  Y_std[miss] <- get("Y_mis")
  resid <- as.vector(Y_std) - rep(get("beta_scaled"), times = rep(n, k))

  Sigma <- kronecker(G_std, A) + kronecker(E_std, diag(n))
  L <- chol(Sigma)
  z <- forwardsolve(t(L), resid)
  joint <- sum(stats::dnorm(z, log = TRUE)) - sum(log(diag(L)))
  expect_equal(sum(get("log_lik")), joint, tolerance = 1e-5)
})

test_that("univariate fits take missing records", {
  skip_if_no_cmdstan()
  sim <- simulate_animal(n_fam = 10, size = 6, G = matrix(0.8), E = matrix(0.8),
                         seed = 15)
  y <- sim$Y[, 1]
  y[c(3, 17)] <- NA
  fit <- lmm_animal(y, sim$X, sim$A, iter = 600, chains = 2, cores = 2, seed = 1,
                    refresh = 0, show_messages = FALSE)
  h2 <- as.numeric(fit$draws("h2", format = "draws_matrix"))
  expect_true(all(is.finite(h2)))
  expect_true("Y_mis" %in% fit$metadata()$stan_variables)
})

test_that("factor structures take missing records too", {
  skip_if_no_cmdstan()
  # 10 families of 6 leaves the diagonal-residual model too weakly identified
  # once twenty gaps are added; this size keeps every configuration sampling
  sim <- simulate_animal(n_fam = 20, size = 8, G = G3, E = E3, seed = 16)
  Y <- sim$Y
  set.seed(17)
  for (i in 1:40) {
    j <- sample(ncol(Y), 1)
    if (sum(!is.na(Y[i, ])) > 1) Y[i, j] <- NA
  }
  for (structure in list(
    list(genetic_structure = "factor", residual_structure = "unstructured",
         n_factors = 2),
    list(genetic_structure = "factor", residual_structure = "shared_factor",
         n_factors = 2),
    list(genetic_structure = "unstructured", residual_structure = "diagonal")
  )) {
    fit <- suppressWarnings(lmm_animal(
      Y, sim$X, sim$A, column_decay = 0.5, iter = 600,
      chains = 2, cores = 2, seed = 1, refresh = 0, show_messages = FALSE,
      genetic_structure = structure$genetic_structure,
      residual_structure = structure$residual_structure,
      n_factors = structure$n_factors))
    expect_true(all(is.finite(draws_mean(fit, "G"))))
  }

  fit_p <- lm_phen_cov(Y, sim$X, covariance_structure = "factor", n_factors = 2,
                       iter = 600, chains = 2, cores = 2, seed = 1, refresh = 0,
                       show_messages = FALSE)
  expect_true(all(is.finite(draws_mean(fit_p, "P"))))
})
