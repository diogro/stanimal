# Input checks shared across the package.
#
# Every check follows the same shape: one helper per kind of value, the argument
# name in backticks, and `call. = FALSE` so the message names the argument the
# caller passed rather than the internal call that raised it. Public functions
# validate at the top, before any work, and reuse these helpers instead of
# repeating conditions inline.

`%||%` <- function(a, b) if (is.null(a)) b else a

# A matrix of finite numbers, with NA allowed only where the data can have
# missing values (marker genotypes, incomplete trait records). Vectors become
# one-column matrices.
.check_numeric_matrix <- function(x, arg, allow_na = FALSE) {
  x <- as.matrix(x)
  if (!is.numeric(x) || any(is.infinite(x)) || (!allow_na && anyNA(x)))
    stop("`", arg, "` must contain only finite numeric values",
         if (allow_na) ", or NA for a missing value", ".", call. = FALSE)
  x
}

# Missing-response bookkeeping shared by the model entry points.
#
# Rows are grouped by which traits they observe: a *pattern*. The phenotypic
# models evaluate a pattern's density on the observed subvector of the
# covariance, so they need each pattern's observed columns and a row order that
# makes every pattern contiguous, which lets Stan slice it with a range. The
# animal models cannot do that: with relatedness a missing value is informed by
# every relative's observed traits, so the marginal does not factorize over rows.
# They fill the gaps with parameters instead, and only need the position of each.
#
# `Y` comes back with the gaps set to zero, so that no NA reaches Stan. The
# observed entries are the ones the standardisation and the likelihood use, and
# both skip the rest explicitly rather than relying on zero standing for a value.
#
# allow_empty_rows is for the models that fill the gaps with parameters. Such a
# model can take a record with nothing observed, because its values are still
# informed by the individual's relatives; one that evaluates each pattern's own
# subvector cannot, because an empty pattern has no likelihood to contribute.
.missing_response <- function(Y, allow_empty_rows = FALSE) {
  miss <- is.na(Y)
  n <- nrow(Y)
  k <- ncol(Y)

  no_data_trait <- which(colSums(!miss) == 0L)
  if (length(no_data_trait))
    stop("every trait in `Y` needs at least one observed value; there is none "
         , "for "
         , paste(.trait_labels(Y)[no_data_trait], collapse = ", ")
         , ". A trait with no data cannot be standardized and has nothing to "
         , "estimate.", call. = FALSE)
  no_data_row <- which(rowSums(!miss) == 0L)
  if (length(no_data_row) && !allow_empty_rows)
    stop("every row of `Y` needs at least one observed trait; row(s) "
         , paste(utils::head(no_data_row, 5L), collapse = ", ")
         , if (length(no_data_row) > 5L) ", ..."
         , " have none. Drop those rows before fitting."
         , call. = FALSE)

  # One key per row, then the distinct keys are the patterns.
  key <- apply(!miss, 1L, function(observed) paste(as.integer(observed), collapse = ""))
  pattern_key <- unique(key)
  n_pattern <- length(pattern_key)
  pattern_id <- match(key, pattern_key)
  pattern_obs <- unname(t(vapply(
    pattern_key,
    function(s) as.integer(strsplit(s, "", fixed = TRUE)[[1L]]),
    integer(k))))

  # Rows grouped by pattern: pattern p occupies rows first_row[p] to
  # first_row[p] + n_rows[p] - 1 of the reordered response.
  row_order <- order(pattern_id)
  sorted <- pattern_id[row_order]
  first_row <- which(!duplicated(sorted))
  n_rows <- tabulate(sorted, nbins = n_pattern)

  Y_observed <- Y
  Y_observed[miss] <- 0
  n_missing <- sum(miss)
  mis_pos <- integer(n * k)
  mis_pos[miss] <- seq_len(n_missing)
  dim(mis_pos) <- dim(miss)
  miss_int <- matrix(as.integer(miss), n, k, dimnames = dimnames(miss))

  list(Y = Y_observed,
       miss = miss_int,
       n_missing = n_missing,
       mis_pos = mis_pos,
       mis_row = row(miss)[miss],
       mis_col = col(miss)[miss],
       n_pattern = n_pattern,
       pattern_obs = pattern_obs,
       row_order = row_order,
       n_rows = n_rows,
       first_row = first_row)
}

# Trait names for a message, falling back to positions.
.trait_labels <- function(Y) {
  labels <- colnames(Y) %||% rep("", ncol(Y))
  blank <- is.na(labels) | !nzchar(labels)
  labels[blank] <- paste("column", which(blank))
  labels
}

.check_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x))
    stop("`", arg, "` must be TRUE or FALSE.", call. = FALSE)
  x
}

# A single whole number, at least `min`: counts, numbers of draws, iterations.
.check_whole <- function(x, arg, min = 1) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < min || x != trunc(x))
    stop("`", arg, "` must be one whole number, at least ", min, ".",
         call. = FALSE)
  # beyond this the coercion to integer is NA, which callers would use as a size
  if (x > .Machine$integer.max)
    stop("`", arg, "` must be at most ", .Machine$integer.max, ".", call. = FALSE)
  as.integer(x)
}

.check_count <- function(x, arg) .check_whole(x, arg, 1)

# Initial number of factors when the user asks for a factor model but leaves
# `n_factors` NULL. A sublinear default avoids the quadratic loadings of
# `n_traits / 2` while still allowing more factors as the trait count grows.
# It only sees the number of traits, so it cannot know whether the data support
# that many: `suggest_n_factors()` is the data-driven answer.
.default_n_factors <- function(n_traits) {
  as.integer(min(n_traits - 1L, max(2L, floor(sqrt(n_traits))), 20L))
}

# Shapes of a beta prior: two finite positive numbers.
.check_beta_prior <- function(x, arg = "beta_prior") {
  if (!is.numeric(x) || length(x) != 2L || any(!is.finite(x)) || any(x <= 0))
    stop("`", arg, "` must be two positive numbers, the shapes of a beta prior.",
         call. = FALSE)
  as.numeric(x)
}

# Warn when a caller explicitly supplies an argument that the selected model
# does not use. Defaults are not warnings, because the caller did not choose
# them; the warning is only for a request that would otherwise be silent.
.warn_ignored <- function(supplied, used, arg) {
  if (supplied && !used)
    warning("`", arg, "` is ignored by this model configuration.",
            call. = FALSE)
  invisible(NULL)
}

# A single probability, strictly inside (0, 1).
.check_prob <- function(x, arg) .check_in_range(x, arg, 0, 1, TRUE, TRUE)

# A single finite number above zero: scale parameters such as the LKJ shape.
.check_positive <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0)
    stop("`", arg, "` must be one finite number above 0.", call. = FALSE)
  as.numeric(x)
}

# A single number inside a range, with either end optionally excluded.
.check_in_range <- function(x, arg, low, high, open_low = FALSE,
                            open_high = FALSE) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x)
  if (ok) ok <- if (open_low) x > low else x >= low
  if (ok) ok <- if (open_high) x < high else x <= high
  if (!ok)
    stop("`", arg, "` must be one number in ", if (open_low) "(" else "[", low,
         ", ", high, if (open_high) ")" else "]", ".", call. = FALSE)
  as.numeric(x)
}

# Iterations of a sampler: whole numbers, at least one draw retained.
.check_iterations <- function(iter, warmup, thin = NULL) {
  iter <- .check_whole(iter, "iter", 1)
  warmup <- .check_whole(warmup, "warmup", 0)
  if (warmup >= iter)
    stop("`warmup` must be smaller than `iter` (", iter, ").", call. = FALSE)
  if (!is.null(thin)) {
    thin <- .check_whole(thin, "thin", 1)
    if (warmup + thin > iter)
      stop("`iter` must leave at least one draw after `warmup`: ",
           warmup + thin, " is above ", iter, ".", call. = FALSE)
  }
  invisible(list(iter = iter, warmup = warmup, thin = thin))
}

.check_sampler_size <- function(chains, cores) {
  chains <- .check_count(chains, "chains")
  cores <- .check_count(cores, "cores")
  if (cores > chains)
    stop("`cores` must not exceed `chains`: there is one core per chain at ",
         "most, and `cores` is ", cores, " for ", chains, " chains.",
         call. = FALSE)
  invisible(list(chains = chains, cores = cores))
}

# The eigendecomposition of A that the animal models take as data. The models
# used to compute it themselves in `transformed data`, once per chain; doing it
# here costs one decomposition per fit, and lets a caller that fits many models
# on the same A (`eigen_A`) pay for it once. Eigenvalues come in ascending
# order, Stan's convention, so the rotated records are in the order the models
# used when they decomposed A themselves; eigenvector signs do not matter to a
# Gaussian likelihood.
.relationship_eigen <- function(A) {
  eig <- eigen(A, symmetric = TRUE)
  ascending <- rev(seq_len(nrow(A)))
  list(Q = eig$vectors[, ascending, drop = FALSE],
       lambda = eig$values[ascending])
}

# Validate a precomputed eigendecomposition of the relationship matrix, as
# passed through `eigen_A`. A silently wrong Q or lambda would give a
# plausible-looking fit of the wrong model, so the pair is checked, not
# assumed: it must reproduce A and Q must be orthonormal (a rescaled Q with
# rescaled eigenvalues reproduces A too, but rotates the data by a
# non-orthogonal map). Both are checked on two fixed probe vectors, which
# costs O(N^2) where rebuilding A would cost O(N^3), and uses no random
# numbers.
.check_eigen_A <- function(eigen_A, A) {
  if (!is.list(eigen_A) || is.null(eigen_A$Q) || is.null(eigen_A$lambda))
    stop("`eigen_A` must be a list with elements `Q` and `lambda`.",
         call. = FALSE)
  Q <- as.matrix(eigen_A$Q)
  lambda <- as.numeric(eigen_A$lambda)
  n <- nrow(A)
  if (!identical(dim(Q), c(n, n)))
    stop("`eigen_A$Q` must be a ", n, " by ", n, " matrix, matching `A`.",
         call. = FALSE)
  if (length(lambda) != n)
    stop("`eigen_A$lambda` must have ", n, " elements, matching `A`.",
         call. = FALSE)
  if (anyNA(Q) || anyNA(lambda) || any(!is.finite(Q)) || any(!is.finite(lambda)))
    stop("`eigen_A` contains missing or non-finite values.", call. = FALSE)
  if (min(lambda) < -.eigen_tolerance(lambda))
    stop("`eigen_A$lambda` must be non-negative (floor the eigenvalues at 0).",
         call. = FALSE)

  index <- seq_len(n)
  probes <- cbind(sin(index), cos(1.3 * index))
  Qt_probes <- crossprod(Q, probes)
  orthonormal <- max(abs(Q %*% Qt_probes - probes)) /
    max(1, max(abs(probes)))
  if (orthonormal > 1e-7)
    stop("`eigen_A$Q` is not orthonormal; pass the eigenvectors of `A` ",
         "unscaled, as returned by `eigen()`.", call. = FALSE)
  target <- A %*% probes
  deviation <- max(abs(Q %*% (lambda * Qt_probes) - target)) /
    max(1, max(abs(target)))
  if (deviation > 1e-7)
    stop("`eigen_A` does not reconstruct `A` (relative deviation ",
         format(deviation, digits = 3),
         "); check that each eigenvalue is paired with its eigenvector.",
         call. = FALSE)
  list(Q = Q, lambda = lambda)
}

# A validated relationship matrix and its eigendecomposition, ready for the
# animal models: A matched to the rows of Y and symmetrized, Q and lambda from
# `eigen_A` when supplied (its rows permuted with A's) or computed here, and
# `blend_A` applied to both. `need` is the text of the error for an A that is
# not positive definite.
.animal_relationship <- function(A, Y, eigen_A = NULL, blend_A = FALSE,
                                 need = "for Stan") {
  original_ids <- rownames(as.matrix(A))
  A <- .check_relationship_matrix(A, Y)
  A <- .check_symmetric(A)
  if (is.null(eigen_A)) {
    eig <- .relationship_eigen(A)
  } else {
    # .check_relationship_matrix() reorders A to the rows of Y by name, and
    # P A P' has eigenvectors P Q, so the rows of Q follow A's.
    if (!is.null(original_ids) && !is.null(rownames(A)) &&
        is.matrix(eigen_A$Q) && nrow(eigen_A$Q) == length(original_ids))
      eigen_A$Q <- eigen_A$Q[match(rownames(A), original_ids), , drop = FALSE]
    eig <- .check_eigen_A(eigen_A, A)
  }
  if (min(eig$lambda) < -.eigen_tolerance(eig$lambda))
    stop("`A` must be positive semidefinite.", call. = FALSE)
  if (blend_A) {
    # 0.99 A + 0.01 I shares A's eigenvectors
    A <- 0.99 * A + 0.01 * diag(nrow(A))
    eig$lambda <- 0.99 * eig$lambda + 0.01
  }
  if (min(eig$lambda) <= .eigen_tolerance(eig$lambda))
    stop("`A` must be positive definite ", need, "; set `blend_A = TRUE` to ",
         "regularize it.", call. = FALSE)
  list(A = A, Q = eig$Q, lambda = pmax(eig$lambda, 0))
}
