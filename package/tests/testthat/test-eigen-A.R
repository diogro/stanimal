# The animal models take the eigendecomposition of A as data (Q, lambda)
# instead of decomposing A in `transformed data` once per chain. lmm_animal()
# decomposes A once per fit in R, or uses a caller's `eigen_A`, so a caller
# fitting many models on the same A pays the O(N^3) decomposition once.

animal_models <- c("animalModelUni.stan", "animalModelUniH2.stan",
                   "animalModelG.stan", "animalModelGH2.stan",
                   "animalModelGFactor.stan", "animalModelGSharedFactor.stan")

model_source <- function(file) {
  path <- system.file("stan", file, package = "stanimal")
  if (!nzchar(path)) path <- file.path("../../inst/stan", file)
  path
}

# The univariate model as it was written before Q and lambda became data: it
# declares A and decomposes it itself, the way a user's custom model_file
# still might.
old_style_uni_model <- function() {
  code <- paste(readLines(model_source("animalModelUni.stan"), warn = FALSE),
                collapse = "\n")
  data_lines <- paste0(
    "  matrix[N, N] Q; // eigenvectors of the relationship matrix A, as columns\n",
    "  vector<lower=0>[N] lambda; // matching eigenvalues of A, floored at 0\n")
  stopifnot(grepl(data_lines, code, fixed = TRUE))
  code <- sub(data_lines, "  matrix[N, N] A;\n", code, fixed = TRUE)
  code <- sub("transformed data {\n", paste0(
    "transformed data {\n",
    "  matrix[N, N] Q = eigenvectors_sym(A);\n",
    "  vector[N] lambda;\n",
    "  {\n    vector[N] ev = eigenvalues_sym(A);\n",
    "    for (n in 1:N) lambda[n] = fmax(ev[n], 0);\n  }\n"), code, fixed = TRUE)
  file <- tempfile(fileext = ".stan")
  writeLines(code, file)
  file
}

eigen_list <- function(A) {
  e <- eigen(A, symmetric = TRUE)
  list(Q = e$vectors, lambda = e$values)
}

test_that("no animal model decomposes A itself", {
  for (file in animal_models) {
    code <- readLines(model_source(file), warn = FALSE)
    expect_false(any(grepl("eigenvectors_sym\\(A\\)|eigenvalues_sym\\(A\\)", code)),
                 label = file)
    expect_true(any(grepl("^  matrix\\[N, N\\] Q;", code)), label = file)
    expect_true(any(grepl("^  vector<lower=0>\\[N\\] lambda;", code)),
                label = file)
    expect_false(any(grepl("^  matrix\\[N, N\\] A;", code)), label = file)
  }
})

test_that("the decomposition follows A when it is reordered to match Y", {
  sim <- simulate_animal(n_fam = 4, size = 5,
                         G = matrix(c(1, 0.2, 0.2, 0.9), 2), E = diag(2))
  ids <- paste0("id", seq_len(nrow(sim$A)))
  Y <- sim$Y
  rownames(Y) <- ids
  shuffle <- c(20:11, 1:10)
  A <- sim$A[shuffle, shuffle]
  dimnames(A) <- list(ids[shuffle], ids[shuffle])

  rel <- stanimal:::.animal_relationship(A, Y, eigen_A = eigen_list(A))
  expect_identical(rownames(rel$A), ids)
  rebuilt <- rel$Q %*% (rel$lambda * t(rel$Q))
  expect_equal(unname(rebuilt), unname(rel$A), tolerance = 1e-10)
})

test_that("blend_A is applied to a supplied decomposition", {
  sim <- simulate_animal(n_fam = 4, size = 5,
                         G = matrix(c(1, 0.2, 0.2, 0.9), 2), E = diag(2))
  rel <- stanimal:::.animal_relationship(sim$A, sim$Y,
                                         eigen_A = eigen_list(sim$A),
                                         blend_A = TRUE)
  blended <- 0.99 * sim$A + 0.01 * diag(nrow(sim$A))
  expect_equal(rel$A, blended)
  expect_equal(rel$Q %*% (rel$lambda * t(rel$Q)), blended, tolerance = 1e-10)
})

test_that("the default decomposition is in Stan's ascending order", {
  sim <- simulate_animal(n_fam = 4, size = 5,
                         G = matrix(c(1, 0.2, 0.2, 0.9), 2), E = diag(2))
  rel <- stanimal:::.animal_relationship(sim$A, sim$Y)
  expect_false(is.unsorted(rel$lambda))
  expect_equal(rel$Q %*% (rel$lambda * t(rel$Q)), sim$A, tolerance = 1e-10)
})

test_that("eigen_A is validated before anything is fitted", {
  sim <- simulate_animal(n_fam = 4, size = 5,
                         G = matrix(c(1, 0.2, 0.2, 0.9), 2), E = diag(2))
  e <- eigen_list(sim$A)
  n <- nrow(sim$A)

  expect_error(lmm_animal(sim$Y, sim$X, sim$A, eigen_A = e$Q),
               "must be a list with elements")
  expect_error(lmm_animal(sim$Y, sim$X, sim$A,
                          eigen_A = list(Q = e$Q[, -1], lambda = e$lambda)),
               "matching `A`")
  expect_error(lmm_animal(sim$Y, sim$X, sim$A,
                          eigen_A = list(Q = e$Q, lambda = e$lambda[-1])),
               "matching `A`")
  expect_error(lmm_animal(sim$Y, sim$X, sim$A,
                          eigen_A = list(Q = e$Q, lambda = -e$lambda)),
               "non-negative")
  expect_error(lmm_animal(sim$Y, sim$X, sim$A,
                          eigen_A = list(Q = e$Q, lambda = rev(e$lambda))),
               "does not reconstruct")
  # a rescaled Q with rescaled eigenvalues reproduces A but is not a rotation
  expect_error(lmm_animal(sim$Y, sim$X, sim$A,
                          eigen_A = list(Q = 2 * e$Q, lambda = e$lambda / 4)),
               "not orthonormal")
  expect_error(lmm_animal(sim$Y, sim$X, sim$A,
                          eigen_A = list(Q = matrix(NA_real_, n, n),
                                         lambda = e$lambda)),
               "missing or non-finite")
})

test_that("a supplied decomposition equal to the default changes nothing", {
  skip_if_no_cmdstan()
  sim <- simulate_animal(n_fam = 8, size = 5,
                         G = matrix(c(1, 0.3, 0.3, 0.8), 2), E = diag(2))
  default <- stanimal:::.animal_relationship(sim$A, sim$Y)
  args <- list(sim$Y, sim$X, sim$A, chains = 1, iter = 60, cores = 1,
               seed = 7, init = 0.5, refresh = 0, show_messages = FALSE)

  fit <- do.call(lmm_animal, args)
  fit_e <- do.call(lmm_animal, c(args, list(
    eigen_A = list(Q = default$Q, lambda = default$lambda))))

  expect_identical(posterior::as_draws_matrix(fit_e$draws("G")),
                   posterior::as_draws_matrix(fit$draws("G")))
})

test_that("a custom model that still declares A reproduces the built-in one", {
  skip_if_no_cmdstan()
  sim <- simulate_animal(n_fam = 8, size = 5, G = matrix(0.8), E = matrix(1))
  args <- list(sim$Y[, 1], sim$X, sim$A, chains = 1, iter = 100, cores = 1,
               seed = 3, init = 0.5, refresh = 0, show_messages = FALSE)

  builtin <- suppressWarnings(do.call(lmm_animal, args))
  custom <- suppressWarnings(do.call(lmm_animal, c(args, list(
    model_file = old_style_uni_model()))))

  # The old model decomposes A in Stan, the new one takes R's decomposition in
  # ascending order. Signs of eigenvectors do not enter a Gaussian likelihood,
  # so the two share every log density and gradient, and the sampler
  # trajectory with them.
  expect_identical(posterior::as_draws_matrix(custom$draws("h2")),
                   posterior::as_draws_matrix(builtin$draws("h2")))
})

test_that("eigen_A works for every structure and parametrization", {
  skip_if_no_cmdstan()
  sim <- simulate_animal(n_fam = 6, size = 5,
                         G = matrix(c(1, 0.3, 0.2, 0.3, 0.9, 0.1, 0.2, 0.1, 0.8), 3),
                         E = diag(3))
  e <- eigen_list(sim$A)
  quick <- list(chains = 1, iter = 40, cores = 1, seed = 2, init = 0.5,
                refresh = 0, show_messages = FALSE)
  configs <- list(
    list(Y = sim$Y[, 1]),
    list(Y = sim$Y[, 1], parametrization = "h2"),
    list(Y = sim$Y),
    list(Y = sim$Y, parametrization = "h2"),
    list(Y = sim$Y, genetic_structure = "factor", n_factors = 1),
    list(Y = sim$Y, genetic_structure = "factor", n_factors = 1,
         residual_structure = "shared_factor")
  )
  for (config in configs) {
    fit <- suppressWarnings(do.call(lmm_animal, c(
      list(X = sim$X, A = sim$A, eigen_A = e), config, quick)))
    expect_s3_class(fit, "CmdStanMCMC")
  }
})

test_that("a different but valid decomposition fits the same model", {
  skip_if_no_cmdstan()
  sim <- simulate_animal(n_fam = 6, size = 5,
                         G = matrix(c(1, 0.2, 0.2, 0.9), 2), E = diag(2))
  # R's own descending order: a different record order, the same model
  args <- list(sim$Y, sim$X, sim$A, chains = 1, iter = 300, cores = 1,
               seed = 11, init = 0.5, refresh = 0, show_messages = FALSE)
  fit <- do.call(lmm_animal, args)
  fit_desc <- do.call(lmm_animal, c(args, list(eigen_A = eigen_list(sim$A))))

  # different trajectories, so compare posterior means, not draws
  g <- colMeans(posterior::as_draws_matrix(fit$draws("G")))
  g_desc <- colMeans(posterior::as_draws_matrix(fit_desc$draws("G")))
  expect_equal(g_desc, g, tolerance = 0.2)
})
