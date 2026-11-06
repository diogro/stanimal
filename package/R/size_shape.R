# Gaps in the traits, for the reparametrization.
#
# Size and shape are a linear map of the traits: size is their mean and shape is
# the residual after it, whose principal components are the shape axes. So
# neither step needs the traits individually, and a missing value is not a
# problem in itself. What is a problem is that a row with a gap has no mean, and
# that the axes need a covariance of the shape space. Both are supplied by one
# model: the traits are multivariate normal, and the mean and covariance are
# estimated by maximum likelihood from whatever was measured (EM, which is exact
# here because every pattern is a linear constraint). A gap is then filled with
# its conditional mean given the traits of the same individual, which is the
# completed value the rest of this file works with.
#
# The completion is a point estimate, so `size_shape()` reports the conditional
# standard deviation of every score it completes (`scores_sd`). It is one
# completion, not a multiple imputation: fitting the scores treats them as data,
# and the uncertainty about the gaps does not travel into that fit.

# A trait never observed, or a pair of traits never observed together, leaves
# part of the covariance unidentified. Say which, rather than let EM return NaN.
.check_gap_pattern <- function(Y) {
  observed <- !is.na(Y)
  blank <- which(rowSums(observed) == 0L)
  if (length(blank))
    stop("every individual needs at least one observed trait, and ",
         length(blank), " row(s) have none (rows ",
         paste(utils::head(blank, 5L), collapse = ", "),
         if (length(blank) > 5L) ", ..." else "", "). Drop them, or leave them ",
         "out of `Y`.", call. = FALSE)
  traits <- .trait_labels(Y)
  rare <- which(colSums(observed) < 2L)
  if (length(rare))
    stop("every trait needs at least two observed values for its variance to be ",
         "estimated: ", paste(traits[rare], collapse = ", "), ".",
         call. = FALSE)
  together <- crossprod(observed + 0)
  diag(together) <- Inf
  thin <- which(together < 2L, arr.ind = TRUE)
  if (nrow(thin)) {
    pairs <- unique(t(apply(thin, 1L, function(ij)
      sort(c(traits[ij[1L]], traits[ij[2L]])))))
    shown <- utils::head(seq_len(nrow(pairs)), 3L)
    stop("traits that are never measured on the same individual have no ",
         "estimable covariance: ",
         paste(vapply(shown, function(r) paste(pairs[r, ], collapse = " and "),
                      character(1)), collapse = "; "),
         if (nrow(pairs) > 3L) "; and others" else "", ".", call. = FALSE)
  }
  invisible(TRUE)
}

# Maximum likelihood mean and covariance of the traits under a multivariate
# normal, by EM. The estimates divide by n rather than by n - 1; the completion
# below uses the covariances only as ratios and the axes are eigenvectors, so
# the difference never reaches the output.
.normal_em <- function(Y, tol = 1e-10, max_iter = 500) {
  centre <- colMeans(Y, na.rm = TRUE)
  scale <- apply(Y, 2L, stats::sd, na.rm = TRUE)
  if (any(!is.finite(scale) | scale <= 0))
    stop("every trait must vary among its observed values for EM completion.",
         call. = FALSE)
  Y <- sweep(sweep(Y, 2L, centre, "-"), 2L, scale, "/")
  missing <- is.na(Y)
  n <- nrow(Y)
  k <- ncol(Y)
  key <- apply(missing, 1L, function(row) paste(as.integer(row), collapse = ""))
  patterns <- lapply(unique(key), function(pattern) {
    rows <- which(key == pattern)
    list(rows = rows, missing = which(missing[rows[1L], ]))
  })

  mean <- colMeans(Y, na.rm = TRUE)
  covariance <- diag(k)
  for (iteration in seq_len(max_iter)) {
    total <- numeric(k)
    product <- matrix(0, k, k)
    for (pattern in patterns) {
      gaps <- pattern$missing
      observed <- setdiff(seq_len(k), gaps)
      rows <- pattern$rows
      completed <- matrix(0, length(rows), k)
      completed[, observed] <- Y[rows, observed, drop = FALSE]
      extra <- matrix(0, k, k)
      if (length(gaps)) {
        # E[y_missing | y_observed] and the covariance left over after it
        gain <- solve(covariance[observed, observed, drop = FALSE],
                      covariance[observed, gaps, drop = FALSE])
        centred <- Y[rows, observed, drop = FALSE] -
          matrix(mean[observed], length(rows), length(observed), byrow = TRUE)
        completed[, gaps] <- sweep(centred %*% gain, 2L, mean[gaps], "+")
        extra[gaps, gaps] <- covariance[gaps, gaps, drop = FALSE] -
          covariance[gaps, observed, drop = FALSE] %*% gain
      }
      total <- total + colSums(completed)
      product <- product + crossprod(completed) + length(rows) * extra
    }
    new_mean <- total / n
    new_covariance <- product / n - tcrossprod(new_mean)
    new_covariance <- (new_covariance + t(new_covariance)) / 2
    if (any(!is.finite(new_mean)) || any(!is.finite(new_covariance)))
      stop("EM completion produced non-finite estimates.", call. = FALSE)
    eigenvalues <- eigen(new_covariance, symmetric = TRUE, only.values = TRUE)$values
    if (min(eigenvalues) <= max(eigenvalues) * k * .Machine$double.eps)
      stop("EM completion produced a singular or non-positive-definite ",
           "covariance; check for collinear traits or insufficient data.",
           call. = FALSE)
    converged <- max(abs(new_covariance - covariance),
                     abs(new_mean - mean)) <
      tol * max(1, max(abs(new_covariance)), max(abs(new_mean)))
    mean <- new_mean
    covariance <- new_covariance
    if (converged) {
      result_mean <- centre + scale * mean
      result_covariance <- covariance * outer(scale, scale)
      if (any(!is.finite(result_mean)) || any(!is.finite(result_covariance)))
        stop("EM completion produced non-finite estimates on the response scale.",
             call. = FALSE)
      return(list(mean = result_mean, covariance = result_covariance,
                  iterations = iteration))
    }
  }
  stop("EM completion did not converge after ", max_iter,
       " iterations; check missingness, collinearity and sample size.",
       call. = FALSE)
}

# Each row completed with its conditional mean, and the conditional covariance
# of its gaps, which is what the standard deviation of a completed score needs.
# Observed entries are left exactly as they were.
.conditional_gaps <- function(Y, mean, covariance) {
  missing <- is.na(Y)
  n <- nrow(Y)
  k <- ncol(Y)
  completed <- Y
  residual_covariance <- array(0, c(n, k, k))
  scale <- sqrt(diag(covariance))
  correlation <- covariance / outer(scale, scale)
  for (i in which(rowSums(missing) > 0L)) {
    gaps <- which(missing[i, ])
    observed <- setdiff(seq_len(k), gaps)
    gain <- solve(correlation[observed, observed, drop = FALSE],
                  correlation[observed, gaps, drop = FALSE])
    completed[i, gaps] <- mean[gaps] +
      scale[gaps] * drop(((Y[i, observed] - mean[observed]) /
                            scale[observed]) %*% gain)
    conditional <-
      correlation[gaps, gaps, drop = FALSE] -
      correlation[gaps, observed, drop = FALSE] %*% gain
    conditional <- (conditional + t(conditional)) / 2
    eigenvalues <- eigen(conditional, symmetric = TRUE, only.values = TRUE)$values
    tolerance <- k * .Machine$double.eps * 10
    if (any(!is.finite(eigenvalues)) || min(eigenvalues) < -tolerance)
      stop("EM completion produced an invalid conditional covariance for row ",
           i, ".", call. = FALSE)
    residual_covariance[i, gaps, gaps] <- conditional * outer(scale[gaps], scale[gaps])
  }
  list(Y = completed, covariance = residual_covariance)
}

#' Share of variance explained by the first principal component
#'
#' Diagnostic for the size and shape reparametrization: when a single component
#' dominates, the traits are mostly one size gradient plus noise, and fitting
#' size and shape scores is a much better behaved model than fitting all of the
#' raw traits.
#'
#' @param Y Numeric matrix of traits, individuals in rows. \code{NA} is allowed
#'   and is not discarded: the correlations are estimated from the available
#'   values, as described in \code{\link{size_shape}}.
#' @return The share of the total variance taken by the first principal
#'   component of the correlation matrix of \code{Y}, between \code{1 / ncol(Y)}
#'   and 1 (1 for a single trait).
#' @details The correlation matrix is used so that traits measured on different
#'   scales do not dominate through their units alone. A share near one means
#'   the traits are one gradient, which is what \code{\link{size_shape}} is for;
#'   a share close to \code{1 / ncol(Y)} means the traits carry independent
#'   information and no reparametrization is needed.
#' @seealso \code{\link{size_shape}}, \code{\link{animal_design}}
#' @examples
#' set.seed(1)
#' G <- random_covariance(3, method = "factor", integration = 0.8)
#' Y <- matrix(rnorm(200 * 3), 200, 3) %*% chol(G)
#' pc1_variance_share(Y)                     # high: one dominant size axis
#' pc1_variance_share(matrix(rnorm(600), 200, 3))  # close to 1 / 3
#' @export
pc1_variance_share <- function(Y) {
  Y <- as.matrix(Y)
  if (!is.numeric(Y))
    stop("`Y` must be a numeric matrix of traits.", call. = FALSE)
  if (any(is.infinite(Y)))
    stop("`Y` must not contain infinite values; use `NA` for a trait that was ",
         "not measured.", call. = FALSE)
  k <- ncol(Y)
  if (is.null(k) || k < 1L)
    stop("`Y` must have at least one column of traits.", call. = FALSE)
  if (k == 1L)
    return(1)
  if (nrow(Y) < 3L)
    stop("`Y` must have at least three individuals.", call. = FALSE)
  if (anyNA(Y)) {
    .check_gap_pattern(Y)
    correlation <- stats::cov2cor(.normal_em(Y)$covariance)
  } else {
    if (any(apply(Y, 2L, stats::sd) == 0))
      stop("every trait in `Y` must vary.", call. = FALSE)
    correlation <- stats::cor(Y)
  }
  eigenvalues <- eigen(correlation, symmetric = TRUE, only.values = TRUE)$values
  eigenvalues[1L] / sum(eigenvalues)
}

#' Size and shape scores for many correlated traits
#'
#' Rewrites a set of correlated traits as log size plus orthogonal shape axes,
#' the reparametrization used to fit many allometric traits: raw traits are
#' usually dominated by a single size gradient, so their genetic covariance is a
#' long, thin ridge that is slow to sample. Size and shape span the same space,
#' so nothing is lost, but the sampler works on well-separated axes.
#'
#' @param Y Numeric matrix of traits, individuals in rows and traits in columns.
#'   Log-transformed measurements are the usual input. Missing values
#'   (\code{NA}) are allowed and are not discarded: see Details.
#' @param n_shape Number of shape axes to keep, between 1 and
#'   \code{ncol(Y) - 1}, and no more than the numerical rank of the centered
#'   shape data. The default requests all \code{ncol(Y) - 1} axes and errors
#'   if the data cannot identify them; choose fewer explicitly in that case.
#'   Keeping all axes makes the transformation a bijection; fewer axes give a
#'   lower-rank (singular) covariance in the original trait space.
#' @return A list of class \code{size_shape} with
#'   \item{scores}{matrix with the size and shape scores, one column each, named
#'     \code{size}, \code{shape1}, \code{shape2}, ...;}
#'   \item{loadings}{\code{ncol(Y)} by \code{1 + n_shape} matrix \code{B} such
#'     that, up to a constant absorbed by the model intercepts,
#'     \code{Y = scores \%*\% t(B)};}
#'   \item{traits}{names of the model traits (\code{size}, \code{shape1}, ...);}
#'   \item{raw_traits}{names of the original traits;}
#'   \item{centre}{column means removed from the shape axes;}
#'   \item{gaps}{number of entries that were not measured, and so were
#'     completed before the transformation (\code{0} for a complete response);}
#'   \item{scores_sd}{standard deviation of each completed score, on the same
#'     layout as \code{scores}, or \code{NULL} when nothing was completed.
#'     Observed individuals have \code{0}.}
#' @details Size is the mean of the traits, \code{rowMeans(Y)}, i.e. the log
#'   geometric mean for log measurements. Shape is \code{Y - size}, which sums
#'   to zero within each individual and therefore lives in an
#'   \code{ncol(Y) - 1} dimensional subspace; its principal components are the
#'   shape axes. Size is often close to the leading axis of variation of a
#'   morphometric data set, but it need not be: it is the unweighted mean of the
#'   traits, while the principal components are shaped by their covariance
#'   structure. The shape axes are ordered by the residual variance they
#'   explain. Their numerical rank is at most \code{min(nrow(Y) - 1,
#'   ncol(Y) - 1)} and can be smaller if the shape data are collinear. It is
#'   computed after completing missing values, using a relative singular-value
#'   tolerance of \code{max(dim(Y)) * .Machine$double.eps}.
#'
#'   The transformation is linear, so genetic and residual covariance matrices
#'   transform as \eqn{B M B^T}; use \code{\link{back_transform}} to do this for
#'   the draws of a fit. Signs of the shape axes are arbitrary.
#'
#'   The prior is not invariant under this change of traits, and that is usually
#'   an asset. The LKJ prior that \code{\link{lmm_animal}} puts on the
#'   correlation matrix acts on correlations between size and shape instead of
#'   correlations between raw traits. The shape axes are principal components of
#'   \code{Y - size}, so they are uncorrelated in the sample and their genetic
#'   correlations are expected to be small; a larger \code{lkj_prior} than for
#'   the same traits fitted raw (stronger shrinkage towards zero correlation) is
#'   therefore reasonable, and it can also make the model cheaper to sample.
#'
#'   \strong{Missing values.} Size is a mean of all the traits and shape is a
#'   difference between them, so a row with a gap has neither. The gap is
#'   therefore filled first: the mean and covariance of the traits are
#'   estimated by maximum likelihood from the values that were measured (an EM
#'   fit of a multivariate normal, which needs no other information), and a
#'   missing entry is replaced by its conditional mean given the traits that
#'   individual does have. The rest of the transformation then runs unchanged,
#'   which is exact for the observed values and a conditional estimate for the
#'   completed ones.
#'
#'   Because those completed scores are estimates, \code{scores_sd} reports the
#'   standard deviation each one carries, and it is worth reading before
#'   comparing a mouse that was missing a trait with one that was not: an entry
#'   completed from a single measured trait is mostly the model's prediction. A
#'   row with nothing observed carries no information and is an error, and two
#'   traits never measured on the same individual make their covariance
#'   unidentifiable, which is also an error.
#'   The EM fit starts from a positive-definite diagonal covariance and errors
#'   if it fails to converge or produces a singular covariance. Check the
#'   missingness pattern, collinearity and sample size in that case.
#'
#'   The completion is a point estimate, not a multiple imputation: the scores
#'   are then treated as data by the model they are passed to, and the
#'   uncertainty in \code{scores_sd} does not travel into that fit. It is the
#'   transformation that is estimated once, exactly as it is for a complete
#'   response.
#' @seealso \code{\link{back_transform}}, \code{\link{animal_design}},
#'   \code{\link{imputed_data}}
#' @examples
#' set.seed(1)
#' G <- random_covariance(3, method = "factor", integration = 0.8)
#' Y <- matrix(rnorm(30 * 3), 30, 3) %*% chol(G)
#' colnames(Y) <- c("a", "b", "c")
#' ss <- size_shape(Y)
#' head(ss$scores)
#' ss$loadings
#'
#' # a gap in one trait is completed from the other two
#' Y[2, "b"] <- NA
#' size_shape(Y)$scores[2, ]
#' @export
size_shape <- function(Y, n_shape = NULL) {
  Y <- as.matrix(Y)
  if (!is.numeric(Y))
    stop("`Y` must be a numeric matrix of traits.", call. = FALSE)
  if (any(is.infinite(Y)))
    stop("`Y` must not contain infinite values; use `NA` for a trait that was ",
         "not measured.", call. = FALSE)
  k <- ncol(Y)
  if (is.null(k) || k < 2L)
    stop("`Y` must have at least two traits.", call. = FALSE)
  if (nrow(Y) < 3L)
    stop("`Y` must have at least three individuals.", call. = FALSE)
  if (is.null(n_shape))
    n_shape <- k - 1L
  n_shape <- .check_count(n_shape, "n_shape")
  if (n_shape > k - 1L)
    stop("`n_shape` must be at most ncol(Y) - 1 (", k - 1L, ").",
         call. = FALSE)

  gaps <- is.na(Y)
  completed <- NULL
  if (any(gaps)) {
    .check_gap_pattern(Y)
    fit <- .normal_em(Y)
    completed <- .conditional_gaps(Y, fit$mean, fit$covariance)
    Y <- completed$Y
  }

  size <- rowMeans(Y)
  shape <- Y - size
  # Shift before centering to avoid spurious variation in identical rows.
  pca <- stats::prcomp(sweep(shape, 2L, shape[1L, ], "-"))
  rank <- min(nrow(Y) - 1L, k - 1L,
              sum(pca$sdev > pca$sdev[1L] *
                    max(dim(shape)) * .Machine$double.eps))
  if (n_shape > rank)
    stop("`n_shape` requests ", n_shape, " axes, but only ", rank,
         " identifiable shape axes are available in the centered data.",
         if (rank > 0L) paste0(" Set `n_shape` to at most ", rank, ".")
         else " The shape data have no varying axis.",
         call. = FALSE)
  axes <- seq_len(n_shape)
  B_shape <- pca$rotation[, axes, drop = FALSE]
  B <- cbind(size = 1, B_shape)

  raw_traits <- colnames(Y) %||% paste0("trait", seq_len(k))
  scores <- cbind(size, pca$x[, axes, drop = FALSE])
  colnames(scores) <- c("size", paste0("shape", axes))
  rownames(B) <- raw_traits
  rownames(scores) <- rownames(Y)

  # size is the mean of the traits and each shape axis is a zero-sum combination
  # of them, so every model trait is a linear function of the k raw traits. That
  # map turns the conditional covariance of a row's gaps into the uncertainty of
  # its scores.
  scores_sd <- NULL
  if (!is.null(completed)) {
    map <- t(cbind(rep(1 / k, k), B_shape))
    scores_sd <- matrix(0, nrow(Y), ncol(scores))
    for (i in which(rowSums(gaps) > 0L)) {
      conditional <- map %*% completed$covariance[i, , ] %*% t(map)
      scores_sd[i, ] <- sqrt(pmax(diag(conditional), 0))
    }
    dimnames(scores_sd) <- dimnames(scores)
  }

  structure(list(scores = scores, loadings = B, traits = colnames(scores),
                 raw_traits = raw_traits, centre = colMeans(shape),
                 gaps = if (any(gaps)) sum(gaps) else 0L,
                 scores_sd = scores_sd),
            class = "size_shape")
}

#' @export
print.size_shape <- function(x, ...) {
  cat("<size_shape>", nrow(x$scores), "individuals,",
      length(x$traits), "model traits from", length(x$raw_traits),
      "raw traits\n")
  cat("model traits:", paste(x$traits, collapse = ", "), "\n")
  cat("raw traits:  ", paste(x$raw_traits, collapse = ", "), "\n")
  if (is.null(x$scores_sd))
    cat("the response was complete\n")
  else
    cat(x$gaps, " entries were not measured and were completed from the trait ",
        "covariance;\n", sum(rowSums(x$scores_sd > 0) > 0),
        " individual(s) have a completed score, with the spread in $scores_sd\n",
        sep = "")
  invisible(x)
}

#' Back-transform posterior draws to the original trait scale
#'
#' Takes a fit of the size and shape scores made with
#' \code{\link{lmm_animal}} and returns draws of the genetic and residual
#' covariance matrices, the phenotypic covariance, the genetic correlations and
#' the heritabilities on the scale of the original traits.
#'
#' @param fit A fit returned by \code{\link{lmm_animal}} (a
#'   \code{CmdStanMCMC} object), or a \code{posterior} draws object, or a matrix
#'   of draws, containing the variables \code{G} and \code{E}, or \code{P} for a
#'   phenotypic fit from \code{\link{lm_phen_cov}}.
#' @param transform A \code{size_shape} object from \code{\link{size_shape}}.
#' @return A list of class \code{size_shape_draws} holding one element per
#'   quantity, in the same layout as \code{\link{lmm_animal_gibbs}}: traits
#'   first and draws last, with every dimension named after the traits. For an
#'   animal fit it holds \code{G}, \code{E}, \code{P}, \code{corrG} and
#'   \code{corrP} as \code{traits} by \code{traits} by \code{draws} arrays,
#'   and \code{h2} as a \code{traits} by \code{draws} matrix. For a
#'   phenotypic fit it holds \code{P} and \code{corrP} as arrays and
#'   \code{sigma_P} as a matrix. \code{traits} holds the trait names.
#' @details For each draw the covariance matrices are mapped as
#'   \eqn{M \to B M B^T}, which is exact because size and shape are a linear
#'   reparametrization of the traits. Heritabilities and correlations are then
#'   recomputed from the transformed matrices: unlike the covariance, they are
#'   not invariant under a change of traits, so the h2 of a shape axis is not
#'   the h2 of any original trait.
#'
#'   The fixed effects \code{beta} stay on the model scale, since a coefficient
#'   for size or a shape axis does not correspond to a single original trait.
#'
#'   Draws are the last dimension, as in \code{\link{lmm_animal_gibbs}}, so
#'   posterior summaries are direct: \code{draws_mean(post, "G")} is the
#'   posterior mean of \code{G} on the trait scale, and
#'   \code{draws_quantile(post, "corrG", probs = 0.05)} its 5\% quantile.
#'   \code{\link{fit_draws}} also accepts this object and returns the same
#'   arrays.
#' @seealso \code{\link{size_shape}}, \code{\link{rename_beta}},
#'   \code{\link{to_draws_list}}
#' @examples
#' ss <- size_shape(cbind(a = rnorm(20), b = rnorm(20), c = rnorm(20)))
#' k <- length(ss$raw_traits)
#' names_mat <- function(var) sprintf(paste0(var, "[%d,%d]"),
#'                                    rep(seq_len(k), k), rep(seq_len(k), each = k))
#' # a real fit supplies G and E; here the same draw is repeated
#' G <- random_covariance(k, min_var = 0.5, max_var = 1.5)
#' E <- random_covariance(k, min_var = 0.5, max_var = 1.5)
#' m <- cbind(matrix(rep(as.vector(G), 4), 4, byrow = TRUE),
#'            matrix(rep(as.vector(E), 4), 4, byrow = TRUE))
#' colnames(m) <- c(names_mat("G"), names_mat("E"))
#' post <- back_transform(m, ss)
#' dim(post$G)                    # traits x traits x draws
#' draws_mean(post, "G")          # posterior mean G on the trait scale
#' draws_mean(post, "h2")
#' @export
back_transform <- function(fit, transform) {
  if (!inherits(transform, "size_shape"))
    stop("`transform` must be a `size_shape` object from `size_shape()`.",
         call. = FALSE)
  B <- transform$loadings
  k_model <- ncol(B)
  k <- nrow(B)

  draws <- if (inherits(fit, "CmdStanFit") || (is.list(fit) && is.function(fit$draws))) {
    variables <- if (is.function(fit$metadata)) fit$metadata()$stan_variables else NULL
    wanted <- intersect(c("G", "E", "P"), variables %||% c("G", "E", "P"))
    if (length(wanted)) fit$draws(wanted) else fit$draws()
  } else fit
  m <- if (inherits(draws, "draws")) unclass(posterior::as_draws_matrix(draws))
       else as.matrix(draws)
  if (!is.matrix(m) || is.null(colnames(m)))
    stop("`fit` must be a fit, a posterior draws object or a named matrix of draws.",
         call. = FALSE)

  in_g <- .matrix_variable_names("G", k_model)
  in_e <- .matrix_variable_names("E", k_model)
  in_p <- .matrix_variable_names("P", k_model)
  has_components <- all(in_g %in% colnames(m)) && all(in_e %in% colnames(m))
  if (!has_components && !all(in_p %in% colnames(m)))
    stop("the draws must contain G and E (lmm_animal) or P (lm_phen_cov) for ",
         k_model, " model traits, as named by the model functions for the ",
         "size_shape scores.", call. = FALSE)

  traits <- transform$raw_traits
  if (has_components) {
    components <- list(G = .transform_draws(m, "G", B, traits),
                       E = .transform_draws(m, "E", B, traits))
    components$P <- components$G + components$E
    components$corrG <- .corr_draws(components$G)
    components$corrP <- .corr_draws(components$P)
    components$h2 <- .diag_draws(components$G) / .diag_draws(components$P)
  } else {
    components <- list(P = .transform_draws(m, "P", B, traits))
    components$corrP <- .corr_draws(components$P)
    components$sigma_P <- sqrt(.diag_draws(components$P))
  }
  structure(c(components, list(traits = traits)), class = "size_shape_draws")
}

#' Convert the back-transformed draws to posterior draws
#'
#' Reshapes the arrays returned by \code{\link{back_transform}} into
#' \code{\link[posterior]{draws_matrix}} objects, one per quantity, so that
#' the \pkg{posterior} and \pkg{bayesplot} tools can be used on them.
#'
#' @param x A \code{size_shape_draws} object from
#'   \code{\link{back_transform}}.
#' @return A named list of \code{draws_matrix} objects, one per quantity of
#'   \code{x}, with the variables named as the models report them:
#'   \code{G[i,j]}, \code{E[i,j]} and \code{P[i,j]} for every element of the
#'   covariance matrices, \code{corrG[i,j]} and \code{corrP[i,j]} for the
#'   correlations above the diagonal only, and \code{h2[i]} or
#'   \code{sigma_P[i]} for the vectors. The trait names are kept in
#'   \code{attr(, "traits")}.
#' @details Correlations are reported once per pair, for \code{i < j}, because
#'   the two halves of the matrix are identical and a draws object has no
#'   notion of symmetry.
#'
#'   To get everything back in a single draws object, bind the elements with
#'   \code{do.call(posterior::bind_draws, unname(draws))}. The names of the
#'   list have to be dropped, because \code{bind_draws()} takes the draws
#'   objects as unnamed arguments.
#' @seealso \code{\link{back_transform}}
#' @examples
#' ss <- size_shape(cbind(a = rnorm(20), b = rnorm(20), c = rnorm(20)))
#' k <- length(ss$raw_traits)
#' G <- random_covariance(k, min_var = 0.5, max_var = 1.5)
#' E <- random_covariance(k, min_var = 0.5, max_var = 1.5)
#' m <- cbind(matrix(rep(as.vector(G), 3), 3, byrow = TRUE),
#'            matrix(rep(as.vector(E), 3), 3, byrow = TRUE))
#' colnames(m) <- c(sprintf("G[%d,%d]", rep(seq_len(k), k), rep(seq_len(k), each = k)),
#'                  sprintf("E[%d,%d]", rep(seq_len(k), k), rep(seq_len(k), each = k)))
#' draws <- to_draws_list(back_transform(m, ss))
#' names(draws)
#' posterior::summarise_draws(draws$corrG)
#' combined <- do.call(posterior::bind_draws, unname(draws))
#' posterior::variables(combined)[1:4]
#' @export
to_draws_list <- function(x) {
  if (!inherits(x, "size_shape_draws"))
    stop("`x` must be a `size_shape_draws` object from `back_transform()`.",
         call. = FALSE)

  components <- setdiff(names(x), "traits")
  out <- lapply(components, function(component) {
    values <- x[[component]]
    if (length(dim(values)) == 3L) {
      k <- dim(values)[1]
      n <- dim(values)[3]
      if (startsWith(component, "corr")) {
        pairs <- which(upper.tri(diag(k)), arr.ind = TRUE)
        pairs <- pairs[order(pairs[, 1], pairs[, 2]), , drop = FALSE]
        index <- cbind(rep(pairs[, 1], each = n), rep(pairs[, 2], each = n),
                       rep(seq_len(n), times = nrow(pairs)))
        draws <- matrix(values[index], n, nrow(pairs))
        variables <- sprintf("%s[%d,%d]", component, pairs[, 1], pairs[, 2])
      } else {
        draws <- t(matrix(values, nrow = k * k))
        variables <- sprintf("%s[%d,%d]", component, rep(seq_len(k), k),
                             rep(seq_len(k), each = k))
      }
    } else {
      draws <- t(values)
      variables <- sprintf("%s[%d]", component, seq_along(x$traits))
    }
    colnames(draws) <- variables
    posterior::as_draws_matrix(draws)
  })
  names(out) <- components
  attr(out, "traits") <- x$traits
  out
}

#' @export
print.size_shape_draws <- function(x, ...) {
  components <- setdiff(names(x), "traits")
  dims <- dim(x[[components[1L]]])
  cat("<size_shape_draws>", dims[length(dims)], "draws,",
      length(x$traits), "traits\n")
  cat("components:", paste(components, collapse = ", "), "\n")
  if (length(x$traits) > 6L) {
    cat("traits:    ", paste(utils::head(x$traits, 6L), collapse = ", "), "...\n")
  } else {
    cat("traits:    ", paste(x$traits, collapse = ", "), "\n")
  }
  invisible(x)
}

# B M B' for every draw of a matrix-valued variable, as trait x trait x draw.
.transform_draws <- function(draws, name, B, traits) {
  k_model <- ncol(B)
  parsed <- array(t(draws[, .matrix_variable_names(name, k_model), drop = FALSE]),
                  c(k_model, k_model, nrow(draws)))
  k <- nrow(B)
  out <- array(apply(parsed, 3L, function(M) B %*% M %*% t(B)), c(k, k, nrow(draws)))
  dimnames(out) <- list(traits, traits, NULL)
  out
}

# cov2cor() of every draw of a trait x trait x draw array.
.corr_draws <- function(A) {
  out <- array(apply(A, 3L, stats::cov2cor), dim(A), dimnames = dimnames(A))
  out
}

# Diagonal of every draw, as a trait x draw matrix.
.diag_draws <- function(A) {
  k <- dim(A)[1]
  n <- dim(A)[3]
  index <- cbind(rep(seq_len(k), times = n), rep(seq_len(k), times = n),
                 rep(seq_len(n), each = k))
  matrix(A[index], k, n, dimnames = list(dimnames(A)[[1]], NULL))
}

# Names of a matrix-valued Stan variable in column-major order, as cmdstanr
# reports them.
.matrix_variable_names <- function(name, k) {
  sprintf("%s[%d,%d]", name, rep(seq_len(k), k), rep(seq_len(k), each = k))
}
