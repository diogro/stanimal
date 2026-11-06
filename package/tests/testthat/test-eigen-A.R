# `eigen_A`: passing the eigendecomposition of A instead of A itself, so that a
# caller fitting many models on the same relationship matrix pays the O(N^3)
# decomposition once instead of once per chain per fit. The variant model is
# animalModelGQR.stan.

eigen_parts <- function(A) {
  eig <- eigen(A, symmetric = TRUE)
  list(Q = eig$vectors[, rev(seq_len(ncol(eig$vectors))), drop = FALSE],
       lambda = pmax(rev(eig$values), 0))
}

model_body <- function(file) {
  lines <- readLines(file, warn = FALSE)
  start <- grep("^model \\{$", lines)[1]
  paste(lines[start:length(lines)], collapse = "\n")
}

test_that("the QR model differs from the original only in how A arrives", {
  plain <- system.file("stan", "animalModelG.stan", package = "stanimal")
  qr <- system.file("stan", "animalModelGQR.stan", package = "stanimal")
  skip_if(!nzchar(plain) || !nzchar(qr), "model sources not available")

  # Everything from the likelihood onwards must be untouched: if it is not, the
  # two models cannot be fitting the same thing.
  expect_identical(model_body(qr), model_body(plain))

  qr_data <- readLines(qr, warn = FALSE)
  start <- grep("^data \\{$", qr_data)[1]
  rest <- qr_data[start:length(qr_data)]
  end <- start + grep("^\\}$", rest)[1] - 1L
  data_block <- qr_data[start:end]
  expect_true(any(grepl("matrix\\[N, N\\] Q;", data_block)))
  expect_true(any(grepl("vector\\[N\\] lambda;", data_block)))
  expect_false(any(grepl("matrix\\[N, N\\] A;", data_block)))
  # and the decomposition is no longer computed in the model
  expect_false(any(grepl("eigenvectors_sym\\(A\\)|eigenvalues_sym\\(A\\)", qr_data)))
})

test_that("eigen_A reproduces the model that takes A", {
  skip_if_no_cmdstan()
  sim <- simulate_animal(n_fam = 8, size = 5,
                         G = matrix(c(1, 0.3, 0.3, 0.8), 2), E = diag(2))
  ep <- eigen_parts(sim$A)

  # Same seed and the same fixed start: two models that share a log density and
  # a gradient share the whole sampler trajectory, so the draws agree bit for
  # bit. A tolerance here would hide a model that is merely close.
  fit_a <- lmm_animal(sim$Y, sim$X, sim$A, chains = 1, iter = 60, cores = 1,
                      seed = 7, init = 0.5, refresh = 0, show_messages = FALSE)
  fit_qr <- lmm_animal(sim$Y, sim$X, sim$A, chains = 1, iter = 60, cores = 1,
                       seed = 7, init = 0.5, eigen_A = ep,
                       refresh = 0, show_messages = FALSE)

  expect_true(attr(fit_qr, "eigen_A"))
  expect_false(attr(fit_a, "eigen_A"))
  expect_identical(posterior::as_draws_matrix(fit_qr$draws("G")),
                   posterior::as_draws_matrix(fit_a$draws("G")))
})

test_that("the eigenvalue order does not matter as long as Q and lambda pair up", {
  skip_if_no_cmdstan()
  sim <- simulate_animal(n_fam = 6, size = 5,
                         G = matrix(c(1, 0.2, 0.2, 0.9), 2), E = diag(2))
  ascending <- eigen_parts(sim$A)
  # R's own (descending) order, unreversed: a different rotation and a different
  # record order, the same model.
  eig <- eigen(sim$A, symmetric = TRUE)
  descending <- list(Q = eig$vectors, lambda = pmax(eig$values, 0))

  fit_asc <- lmm_animal(sim$Y, sim$X, sim$A, chains = 1, iter = 300, cores = 1,
                        seed = 11, init = 0.5, eigen_A = ascending,
                        refresh = 0, show_messages = FALSE)
  fit_desc <- lmm_animal(sim$Y, sim$X, sim$A, chains = 1, iter = 300, cores = 1,
                         seed = 11, init = 0.5, eigen_A = descending,
                         refresh = 0, show_messages = FALSE)

  # A different record order means different trajectories, so compare the
  # posterior means rather than the draws.
  g_asc <- colMeans(posterior::as_draws_matrix(fit_asc$draws("G")))
  g_desc <- colMeans(posterior::as_draws_matrix(fit_desc$draws("G")))
  expect_equal(g_desc, g_asc, tolerance = 0.2)
})

test_that("eigen_A is validated before anything is fitted", {
  sim <- simulate_animal(n_fam = 4, size = 5,
                         G = matrix(c(1, 0.2, 0.2, 0.9), 2), E = diag(2))
  ep <- eigen_parts(sim$A)
  n <- nrow(sim$A)

  expect_error(lmm_animal(sim$Y, sim$X, sim$A, eigen_A = ep$Q),
               "must be a list with elements")
  expect_error(lmm_animal(sim$Y, sim$X, sim$A,
                          eigen_A = list(Q = ep$Q[, -1], lambda = ep$lambda)),
               "matching `A`")
  expect_error(lmm_animal(sim$Y, sim$X, sim$A,
                          eigen_A = list(Q = ep$Q, lambda = ep$lambda[-1])),
               "matching `A`")
  expect_error(lmm_animal(sim$Y, sim$X, sim$A,
                          eigen_A = list(Q = ep$Q, lambda = -ep$lambda)),
               "non-negative")
  expect_error(lmm_animal(sim$Y, sim$X, sim$A,
                          eigen_A = list(Q = ep$Q, lambda = rep(1, n))),
               "does not reconstruct")
  expect_error(lmm_animal(sim$Y, sim$X, sim$A,
                          eigen_A = list(Q = matrix(NA_real_, n, n),
                                         lambda = ep$lambda)),
               "missing or non-finite")
})

test_that("eigen_A is refused where it cannot be honoured", {
  sim <- simulate_animal(n_fam = 4, size = 5,
                         G = matrix(c(1, 0.2, 0.2, 0.9), 2), E = diag(2))
  ep <- eigen_parts(sim$A)

  expect_error(lmm_animal(sim$Y, sim$X, sim$A, eigen_A = ep,
                          parametrization = "h2"),
               "parametrization")
  expect_error(lmm_animal(sim$Y, sim$X, sim$A, eigen_A = ep,
                          genetic_structure = "factor"),
               "factor")
  expect_error(lmm_animal(sim$Y, sim$X, sim$A, eigen_A = ep,
                          model_file = tempfile()),
               "model_file")
  expect_error(lmm_animal(sim$Y[, 1], sim$X, sim$A,
                          eigen_A = list(Q = ep$Q, lambda = ep$lambda)),
               "multivariate")
})
