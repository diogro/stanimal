scree_modules <- function(n = 400, p = 6, seed = 1) {
  loadings <- matrix(0, p, 2)
  loadings[seq_len(p / 2), 1] <- 0.8
  loadings[p / 2 + seq_len(p / 2), 2] <- 0.8
  P <- loadings %*% t(loadings) + diag(p) * 0.36
  set.seed(seed)
  Y <- matrix(stats::rnorm(n * p), n) %*% chol(P)
  colnames(Y) <- paste0("t", seq_len(p))
  Y
}

test_that("scree_eigenvalues matches the eigenvalues of the correlation", {
  Y <- scree_modules()
  scree <- scree_eigenvalues(Y)
  expect_s3_class(scree, "scree")
  expect_equal(names(scree), c("component", "eigenvalue", "share", "cumulative"))
  expect_equal(scree$component, seq_len(ncol(Y)))
  expect_equal(scree$eigenvalue,
               eigen(stats::cor(Y), symmetric = TRUE, only.values = TRUE)$values)
  expect_equal(scree$eigenvalue, sort(scree$eigenvalue, decreasing = TRUE))
  expect_equal(sum(scree$share), 1)
  expect_equal(scree$eigenvalue / ncol(Y), scree$share)
  expect_equal(tail(scree$cumulative, 1), 1)
})

test_that("scree_eigenvalues uses the residuals, not the marginal traits", {
  Y <- scree_modules()
  X <- cbind(Intercept = 1, shift = stats::rnorm(nrow(Y), sd = 3))
  # one strong shared direction from the covariate: the marginal correlation is
  # dominated by it, and removing it must recover the two-module structure
  Y_shifted <- Y + X[, 2] %*% t(rep(1, ncol(Y))) * 2
  colnames(Y_shifted) <- colnames(Y)

  marginal <- scree_eigenvalues(Y_shifted)
  residual <- scree_eigenvalues(Y_shifted, X)
  expect_gt(marginal$eigenvalue[1], residual$eigenvalue[1])
  expect_equal(residual$eigenvalue,
               eigen(stats::cor(stats::lm.fit(X, Y_shifted)$residuals),
                     symmetric = TRUE, only.values = TRUE)$values)
  expect_equal(as.integer(suggest_n_factors(Y_shifted, X)), 2)
  # ignoring the covariate collapses the structure into its one shared direction
  expect_equal(as.integer(suggest_n_factors(Y_shifted)), 1)})

test_that("suggest_n_factors counts the dimensions that clear the noise", {
  Y <- scree_modules()
  set.seed(11)
  expect_equal(as.integer(suggest_n_factors(Y)), 2)

  set.seed(12)
  noise <- matrix(stats::rnorm(400 * 10), 400, 10,
                  dimnames = list(NULL, paste0("z", 1:10)))
  flat <- suggest_n_factors(noise)
  expect_equal(as.integer(flat), 1)           # the models need at least one
  expect_equal(attr(flat, "retained"), 0)     # nothing cleared the reference
  expect_length(attr(flat, "reference"), 10)
  expect_length(attr(flat, "eigenvalues"), 10)
})

test_that("suggest_n_factors respects its caps", {
  Y <- scree_modules()
  set.seed(13)
  expect_lte(suggest_n_factors(Y, max_factors = 1), 1)
  set.seed(14)
  expect_lte(suggest_n_factors(Y), ncol(Y) - 1)
  set.seed(15)
  loose <- suggest_n_factors(Y, nsim = 50, level = 0.9)
  expect_true(loose >= suggest_n_factors(Y))
})

test_that("scree helpers validate their input", {
  Y <- scree_modules()
  expect_error(suggest_n_factors(Y[, 1, drop = FALSE]), "at least two traits")
  expect_error(suggest_n_factors(Y, level = 1), "between 0 and 1")
  expect_error(suggest_n_factors(Y, level = 0), "between 0 and 1")
  expect_error(suggest_n_factors(Y, nsim = 0), "`nsim`")
  expect_error(suggest_n_factors(Y, max_factors = 0), "`max_factors`")

  # a covariate that reproduces a trait leaves no residual variance for it
  X <- cbind(Intercept = 1, copy = Y[, 1])
  expect_error(scree_eigenvalues(Y, X), "no residual variance")
})

test_that("plot.scree draws the eigenvalues and the reference line", {
  skip_if_not_installed("ggplot2")
  Y <- scree_modules()
  scree <- scree_eigenvalues(Y)
  expect_s3_class(plot(scree), "ggplot")
  set.seed(16)
  suggestion <- suggest_n_factors(Y)
  expect_s3_class(plot(scree, reference = suggestion), "ggplot")
  expect_s3_class(plot(scree, reference = attr(suggestion, "reference")), "ggplot")
  expect_error(plot(scree, reference = c(1, 2)), "one value per component")
})
