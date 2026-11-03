# Estimated breeding values from a marginalized fit.
#
# The animal models integrate the additive genetic effects out. That is what
# makes them scale to many individuals and many traits, and it is also why `u`
# is never reported: the sampler works on the marginal likelihood of the
# response, where `u` has already been summed away.
#
# Nothing stops it being recovered afterwards. The model says
#
#   y = X beta + u + e,   u ~ N(0, A (x) G),   e ~ N(0, I (x) E),
#
# and the best linear unbiased predictor of `u` given the data and given the
# covariance parameters is the conditional mean
#
#   u_hat = (A (x) G) V^-1 (y - X beta),   V = A (x) G + I (x) E
#         = Cov(u, y) Var(y)^-1 (y - X beta).
#
# That is a linear function of the response, which is what makes the missing
# data case cheap (see below).
#
# The N x N system is never formed. With A = Q diag(lambda) Q', the rotated
# response Q'y has independent records and record i has covariance
# lambda_i G + E, so
#
#   u_hat_rot[i, ] = lambda_i G (lambda_i G + E)^-1 (Q'y)[i, ],
#
# one k x k solve per individual instead of one N-by-N solve, and rotating back
# is one matrix product. This is the same eigendecomposition the model block
# already uses, recomputed here because it is O(N^3) once and the whole point is
# not to need the sampler's internals.

# The conditional distribution of u given the data, at fixed G, E and beta.
# Returns the predicted values and, optionally, their conditional variance:
#
#   Var(u | y) = A (x) G - (A (x) G) V^-1 (A (x) G),
#
# whose diagonal is what an approximate standard error and a reliability are
# built from. In the rotated basis this is block diagonal, one k x k block per
# record, and the marginal variance of individual n comes back by the same
# rotation: Var(u[n, k]) = sum_i Q[n, i]^2 * Var(u_rot[i, k]).
.conditional_breeding_values <- function(Y, X, A, G, E, beta, se = TRUE) {
  n <- nrow(Y)
  k <- ncol(Y)
  traits <- colnames(Y) %||% paste0("trait", seq_len(k))
  individuals <- rownames(Y) %||% as.character(seq_len(n))

  eig <- eigen(A, symmetric = TRUE)
  # A is positive semidefinite but rounding can leave a tiny negative value
  lambda <- pmax(eig$values, 0)
  Q <- eig$vectors

  residual <- Y - X %*% t(matrix(beta, nrow = k))
  rotated <- crossprod(Q, residual)

  u_rot <- matrix(0, n, k)
  variance <- if (se) matrix(0, n, k) else NULL
  for (i in seq_len(n)) {
    # lambda_i G + E is positive definite whenever E is, which every model
    # guarantees: s_E > 0 times a correlation matrix, or s_E > 0 diagonal.
    precision <- lambda[i] * G + E
    # lambda_i G (lambda_i G + E)^-1, written transposed so that the solve runs
    # with G as the right-hand side rather than inverting first
    gain <- lambda[i] * t(solve(precision, G))
    u_rot[i, ] <- drop(gain %*% rotated[i, ])
    if (se) {
      conditional <- lambda[i] * G - gain %*% precision %*% t(gain)
      variance[i, ] <- diag(conditional)
    }
  }

  values <- Q %*% u_rot
  dimnames(values) <- list(individuals, traits)
  if (!se)
    return(list(values = values, se = NULL, reliability = NULL))

  # The marginal variance of each prediction, rotating the per-record blocks
  marginal <- (Q^2) %*% variance
  standard_error <- sqrt(pmax(marginal, 0))
  dimnames(standard_error) <- list(individuals, traits)
  reliability <- 1 - sweep(marginal, 2L, diag(G), "/")
  # the diagonal of a conditional variance cannot exceed the unconditional one;
  # only rounding can push it past, so clamp rather than report 1 + 1e-16
  reliability <- pmin(pmax(reliability, 0), 1)
  dimnames(reliability) <- list(individuals, traits)
  list(values = values, se = standard_error, reliability = reliability)
}

# The posterior mean of a completed dataset, over the draws imputed_data()
# returns. Averaging is legitimate here, unlike for a downstream analysis on
# imputed data: the predictor below is linear in the response, so the mean of
# the completed datasets gives exactly the conditional mean of u.
.completed_mean <- function(imputation) {
  completed <- imputation$Y
  d <- dim(completed)
  out <- apply(completed, c(1L, 2L), mean)
  matrix(out, nrow = d[1L], dimnames = dimnames(completed)[1:2])
}

#' Estimated breeding values (BLUPs)
#'
#' Predicts the additive genetic value of each individual from a fitted animal
#' model, on the response scale. The model itself integrates the breeding values
#' out and so never reports them; this recovers them by post-processing the fit,
#' which is exact for the conditional mean once \code{G}, \code{E} and
#' \code{beta} are set to their posterior means.
#'
#' @param fit A fit from \code{\link{lmm_animal}}.
#' @param se Logical: also report the conditional standard error of each
#'   prediction and its reliability. Default \code{TRUE}.
#' @param ndraws Number of completed datasets to average over when the response
#'   had missing values; ignored when it is complete.
#' @param seed Optional seed for those completions, so that a run is
#'   reproducible. The sampler's seed does not cover this step.
#' @return An object of class \code{stanimal_breeding_values}: a list with
#'   \code{values} (individuals by traits, on the response scale),
#'   \code{se} and \code{reliability} (the same shape, or \code{NULL} when
#'   \code{se = FALSE}), the \code{G} and \code{E} matrices used, and the
#'   individual and trait labels. \code{print()} shows the predictions of the
#'   first few individuals with their standard errors and reliabilities.
#' @details Write the model as \eqn{y = X\beta + u + e} with
#'   \eqn{u \sim N(0, A \otimes G)} and \eqn{e \sim N(0, I \otimes E)}. The
#'   predicted breeding values are the conditional means
#'   \deqn{\hat u = (A \otimes G) V^{-1} (y - X\beta), \quad V = A \otimes G + I \otimes E,}
#'   evaluated at the posterior means of \code{G}, \code{E} and \code{beta}.
#'   This is therefore an empirical Bayes estimate rather than a posterior: it
#'   does not propagate the uncertainty in the covariance parameters, and the
#'   standard errors below are conditional on them.
#'
#'   \code{se} is the square root of \eqn{\mathrm{Var}(u \mid y)}, which is
#'   \eqn{A \otimes G - (A \otimes G)V^{-1}(A \otimes G)} at the same parameter
#'   values. \code{reliability} is \eqn{1 - \mathrm{Var}(u_{nk} \mid y)/G_{kk}},
#'   the share of a trait's additive variance that the prediction accounts for,
#'   between 0 and 1. A value near 1 means the individual's relatives and its own
#'   records pin its breeding value down; near 0 means the prediction is close to
#'   the population mean. These are the usual BLUP quantities, so they are
#'   directly comparable with what a classical animal-model program reports.
#'
#'   Predictions are deviations from the fixed effects, so they are on the scale
#'   of the trait and can be read as "this individual is above or below the mean
#'   of its group by this much".
#'
#'   Missing values are handled by averaging the completed datasets from
#'   \code{\link{imputed_data}} and predicting from that mean, which is not the
#'   approximation it looks like: the predictor is linear in the response, so
#'   plugging in \eqn{E[y \mid y_\mathrm{observed}]} gives
#'   \eqn{E[u \mid y_\mathrm{observed}]} exactly.
#'
#'   Because \code{A} enters through its eigendecomposition, this is also a
#'   direct check on a fit: with \code{A} equal to the identity the predictions
#'   are ordinary ridge shrinkage of each individual's own residual, and with
#'   relatedness they borrow from relatives.
#' @seealso \code{\link{lmm_animal}}, \code{\link{imputed_data}},
#'   \code{\link{draws_mean}}
#' @examples
#' \dontrun{
#' fit <- lmm_animal(Y, X, A, iter = 1000, warmup = 500,
#'                   chains = 2, cores = 2, refresh = 0)
#' ebv <- breeding_values(fit)
#' ebv                                   # first individuals, with reliability
#' head(ebv$values)                      # individuals x traits
#' range(ebv$reliability)
#'
#' # The predicted values are on the trait's scale, so they can be ranked
#' order(ebv$values[, 1], decreasing = TRUE)[1:5]
#' }
#' @export
#' @author Diogo Melo
breeding_values <- function(fit, se = TRUE, ndraws = 100, seed = NULL) {
  se <- .check_flag(se, "se")
  data <- .stanimal_data(fit, "its breeding values cannot be predicted")
  if (!identical(data$family, "animal"))
    stop("`breeding_values()` needs an animal model: `lm_phen_cov()` has no ",
         "genetic component, so there is no breeding value to predict.",
         call. = FALSE)

  Y <- data$Y
  n <- nrow(Y)
  k <- ncol(Y)
  traits <- .draw_traits(fit, NULL, k)

  # one pass over the draws, then the summaries: a fit that reports G and E
  # reports them on the response scale already, which is the scale the response
  # and the design are on
  d <- fit_draws(fit)
  G <- if ("G" %in% names(d)) draws_mean(d, "G") else NULL
  E <- if ("E" %in% names(d)) draws_mean(d, "E") else NULL
  if (is.null(G) || is.null(E)) {
    # a univariate fit reports the standard deviations, not the covariances
    sigma_G <- if ("sigma_G" %in% names(d)) draws_mean(d, "sigma_G") else NULL
    sigma_E <- if ("sigma_E" %in% names(d)) draws_mean(d, "sigma_E") else NULL
    if (is.null(sigma_G) || is.null(sigma_E))
      stop("this fit reports neither `G` and `E` nor `sigma_G` and `sigma_E`, ",
           "so the breeding values cannot be predicted from it.", call. = FALSE)
    G <- diag(as.numeric(sigma_G)^2, k)
    E <- diag(as.numeric(sigma_E)^2, k)
  }
  G <- matrix(G, k, k, dimnames = list(traits, traits))
  E <- matrix(E, k, k, dimnames = list(traits, traits))

  beta <- if ("beta" %in% names(d)) draws_mean(d, "beta") else NULL
  if (is.null(beta))
    stop("this fit does not report `beta`, so the fixed effects cannot be ",
         "removed before predicting.", call. = FALSE)
  beta <- matrix(beta, nrow = k)
  if (ncol(beta) != ncol(data$X))
    stop("this fit reports ", ncol(beta), " fixed effects but the design it ",
         "was given has ", ncol(data$X), " column(s).", call. = FALSE)

  if (anyNA(Y))
    Y <- .completed_mean(imputed_data(fit, ndraws = ndraws, seed = seed))

  predicted <- .conditional_breeding_values(Y, data$X, data$A, G, E, beta,
                                            se = se)
  structure(list(values = predicted$values, se = predicted$se,
                 reliability = predicted$reliability, G = G, E = E,
                 traits = traits, individuals = rownames(predicted$values)),
            class = "stanimal_breeding_values")
}

#' @rdname breeding_values
#' @param x An object from \code{breeding_values()}.
#' @param n Number of individuals to show. Default 4.
#' @param ... Unused; for method compatibility.
#' @export
print.stanimal_breeding_values <- function(x, n = 4L, ...) {
  values <- x$values
  n <- min(as.integer(n), nrow(values))
  cat("<stanimal_breeding_values>", nrow(values), "individual(s),",
      ncol(values), "trait(s), empirical Bayes\n")
  if (is.null(x$reliability))
    cat("G, E and beta at their posterior means; no standard errors requested\n")
  else
    cat("G, E and beta at their posterior means; reliability is 1 - Var(u | y) / diag(G)\n")

  selected <- seq_len(n)
  out <- data.frame(
    individual = rep(x$individuals[selected], times = ncol(values)),
    trait = rep(colnames(values), each = n),
    value = round(as.numeric(values[selected, , drop = FALSE]), 3),
    check.names = FALSE
  )
  if (!is.null(x$se)) {
    out$se <- round(as.numeric(x$se[selected, , drop = FALSE]), 3)
    out$reliability <- round(as.numeric(x$reliability[selected, , drop = FALSE]), 3)
  }
  print(out, row.names = FALSE)
  invisible(x)
}
