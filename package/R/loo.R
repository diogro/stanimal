# Predictive model comparison for a fitted model.
#
# Phenotypic models report one `log_lik` term per individual, with missing
# traits integrated out, and those go to `loo` directly. Keep the chain
# structure intact for its effective sample size correction.
#
# Animal models are not factorized over individuals: relatives share breeding
# values, so the density of one animal depends on every other. Their `log_lik`
# terms are rotated records (independent after the eigendecomposition of A for
# complete responses), which sum to the joint density but do not each belong to
# one animal, and with missing responses they are completed-data terms.
#
# The default for animal fits is therefore computed here, after the fit, as
# leave-one-individual-out for a non-factorized normal model (Buerkner, Gabry
# and Vehtari 2021, Computational Statistics 36:1243-1261). For one posterior
# draw the term of individual i is
#
#   log p(y_i,obs | y_-i, theta),
#
# the conditional density of i's observed traits given every other individual,
# with i's own gaps integrated out. When other individuals have gaps, the gaps
# are latent parameters of the animal model (`Y_mis`), so the draw conditions on
# their sampled values. Importance sampling with these terms is exact: with phi
# the parameters and the other individuals' gaps, p(phi | y_-i) is proportional
# to p(phi | y) / p(y_i,obs | phi, y_-i,obs), and the right-hand density is the
# term above. Averaging a term over the others' gaps recovers the observed-data
# conditional p(y_i,obs | y_-i,obs, theta).
#
# The conditional comes from the precision of the stacked response. With
# A = Q diag(lambda) Q' and C C' = L_E^-1 G L_E^-T = U diag(d) U', write
# Z = Q' R L_E^-T U for the residuals R (individuals by traits) and
# D[n, k] = 1 + lambda_n d_k. Then, with B = L_E U,
#
#   H = Q (Z / D),   W = (Q * Q) (1 / D),
#   Var(y_i | y_-i) = B diag(1 / W[i, ]) B',
#   y_i - E(y_i | y_-i) = B (H[i, ] / W[i, ]),
#
# so for an individual with every trait observed
#
#   log p = -1/2 (K log 2 pi + 2 sum log diag L_E - sum_k log W[i, k]
#                 + sum_k H[i, k]^2 / W[i, k]),
#
# and one with gaps takes the observed block of the same mean and covariance.
# Per draw that is two N x N by N x K products, batched over draws.
#
# The terms are densities of the standardized response, as the model's own
# `log_lik` is: the response-scale density plus the log standard deviation of
# each observed trait.

.check_loo_data <- function(fit, unit = "record") {
  data <- attr(fit, "stanimal_data")
  if (is.null(data) || !is.list(data) ||
      !(identical(data$family, "animal") ||
        identical(data$family, "phenotypic")) || is.null(data$Y)) {
    warning("Cannot verify whether `log_lik` is suitable for LOO or WAIC: ",
            "response and model metadata are unavailable. Ensure the terms ",
            "are appropriate observed-data predictive units; incomplete animal ",
            "fits' completed-data terms are not.", call. = FALSE)
  } else if (identical(data$family, "animal") && anyNA(data$Y)) {
    stop("rotated-record LOO and WAIC (`unit = \"record\"`) are unavailable ",
         "for animal fits with missing responses: those terms use imputed ",
         "values rather than an observed-data likelihood. Use the default ",
         "`unit = \"individual\"`.", call. = FALSE)
  }
}

# Draws of `log_lik` with their chain structure: a draws by records matrix, as
# `loo` wants, and the chain each draw came from.
.log_lik_draws <- function(fit) {
  if (inherits(fit, "stanimal_laplace"))
    stop("`fit_loo()` needs posterior draws. A `stanimal_laplace` object is a "
         , "Gaussian approximation at the mode, whose draws are not weighted "
         , "posterior draws, so a comparison built from them would be "
         , "meaningless. Fit with `lmm_animal()` or `lm_phen_cov()` first.",
         call. = FALSE)
  if (.is_gibbs_draws(fit))
    stop("`lmm_animal_gibbs()` does not support LOO or WAIC. Use "
         , "`lmm_animal()` for model comparison.", call. = FALSE)
  if (inherits(fit, "stanfit")) {
    .check_loo_data(fit)
    # rstan keeps the declared parameters here, so this asks the fit instead of
    # catching a failed extraction
    draws <- if ("log_lik" %in% methods::slot(fit, "model_pars"))
      as.array(fit, pars = "log_lik") else NULL
  } else if (inherits(fit, "CmdStanFit") ||
             (is.list(fit) && is.function(fit$draws))) {
    .check_loo_data(fit)
    # ask the fit what it holds rather than catching a failure to read it: a
    # genuine read error must not be reported as a missing `log_lik`
    if (!"log_lik" %in% fit$metadata()$stan_variables)
      draws <- NULL
    else
      draws <- as.array(fit$draws(variables = "log_lik",
                                  format = "draws_array"))
  } else if (is.matrix(fit) && !is.null(colnames(fit))) {
    .check_loo_data(fit)
    draws <- fit[, grep("^log_lik", colnames(fit)), drop = FALSE]
    return(list(log_lik = draws, chain_id = rep(1L, nrow(draws))))
  } else {
    stop("`x` must be a fit from `lmm_animal()` or `lm_phen_cov()`, or a draws ",
         "matrix with `log_lik` columns.", call. = FALSE)
  }
  if (is.null(draws) || length(dim(draws)) != 3L || dim(draws)[3L] < 1L)
    stop("this fit does not report `log_lik`, so LOO and WAIC cannot be ",
         "computed. Every built-in model reports it; a custom `model_file` has ",
         "to compute it in its generated quantities.", call. = FALSE)
  list(log_lik = matrix(draws, nrow = dim(draws)[1L] * dim(draws)[2L]),
       chain_id = rep(seq_len(dim(draws)[2L]), each = dim(draws)[1L]))
}

# The draws an individual-level term needs, as a draws matrix with the chain of
# each row. Only the covariance, fixed-effect and gap variables are read: the
# full draws include `log_lik`, one column per record.
.individual_loo_draws <- function(fit) {
  wanted <- c("G", "E", "sigma_G", "sigma_E", "beta", "Y_mis")
  if (inherits(fit, "stanfit")) {
    present <- intersect(wanted, methods::slot(fit, "model_pars"))
    arr <- as.array(fit, pars = present)
  } else {
    present <- intersect(wanted, fit$metadata()$stan_variables)
    arr <- fit$draws(variables = present, format = "draws_array")
  }
  arr <- posterior::as_draws_array(arr)
  n_iter <- posterior::niterations(arr)
  n_chain <- posterior::nchains(arr)
  list(draws = as.matrix(posterior::as_draws_matrix(arr)),
       chain_id = rep(seq_len(n_chain), each = n_iter))
}

# Pointwise log density of each individual's observed traits given every other
# individual, for a batch of draws (see the header of this file). `R` is the
# completed residual, individuals x traits x draws, on the response scale; `G`
# and `E` are traits x traits x draws. Rows with nothing observed get NA.
.individual_log_density <- function(R, G, E, Q, QQ, lambda, miss) {
  n <- dim(R)[1L]
  k <- dim(R)[2L]
  b <- dim(R)[3L]
  rotated <- array(0, c(n, k, b))
  log_det_E <- numeric(b)
  basis <- array(0, c(k, k, b))
  inv_den <- matrix(0, n, k * b)
  for (s in seq_len(b)) {
    E_s <- matrix(E[, , s], k, k)
    G_s <- matrix(G[, , s], k, k)
    L_E <- t(chol(E_s))
    M <- forwardsolve(L_E, t(forwardsolve(L_E, G_s)))
    eig <- eigen((M + t(M)) / 2, symmetric = TRUE)
    d <- pmax(eig$values, 0)
    U <- eig$vectors
    # R L_E^-T U, so that Q' of it is Z
    rotated[, , s] <- t(forwardsolve(L_E, t(matrix(R[, , s], n, k)))) %*% U
    log_det_E[s] <- 2 * sum(log(diag(L_E)))
    basis[, , s] <- L_E %*% U
    inv_den[, (s - 1L) * k + seq_len(k)] <- 1 / (1 + outer(lambda, d))
  }
  Z <- crossprod(Q, matrix(rotated, n, k * b))
  H <- Q %*% (Z * inv_den)
  W <- QQ %*% inv_den
  complete <- rowSums(miss) == 0L
  partial <- which(!complete & rowSums(!miss) > 0L)
  out <- matrix(NA_real_, b, n)
  for (s in seq_len(b)) {
    cols <- (s - 1L) * k + seq_len(k)
    H_s <- H[, cols, drop = FALSE]
    W_s <- W[, cols, drop = FALSE]
    out[s, complete] <- -0.5 * (k * log(2 * pi) + log_det_E[s] -
      rowSums(log(W_s[complete, , drop = FALSE])) +
      rowSums(H_s[complete, , drop = FALSE]^2 / W_s[complete, , drop = FALSE]))
    if (length(partial)) {
      B <- matrix(basis[, , s], k, k)
      for (i in partial) {
        o <- !miss[i, ]
        B_o <- B[o, , drop = FALSE]
        S_o <- B_o %*% (t(B_o) / W_s[i, ])
        e_o <- drop(B_o %*% (H_s[i, ] / W_s[i, ]))
        L <- chol((S_o + t(S_o)) / 2)
        z <- backsolve(L, e_o, transpose = TRUE)
        out[s, i] <- -0.5 * (sum(o) * log(2 * pi) + 2 * sum(log(diag(L))) +
                               sum(z^2))
      }
    }
  }
  out
}

# Leave-one-individual-out terms for an animal fit: a draws by individuals
# matrix of log densities on the standardized response scale, the chain of each
# draw, and the individuals the columns belong to (those with at least one
# observed trait).
.individual_log_lik <- function(fit, eigen_A = NULL, batch = NULL) {
  data <- attr(fit, "stanimal_data")
  if (isTRUE(attr(fit, "custom_model")))
    stop("individual-level LOO assumes the built-in animal likelihood, which a ",
         "custom `model_file` may not use. Use `unit = \"record\"` for the ",
         "model's own `log_lik` terms.", call. = FALSE)
  if (is.null(data) || !identical(data$family, "animal") || is.null(data$A))
    stop("individual-level LOO needs the response, design and relationship ",
         "matrix the fit was given, which fits from `lmm_animal()` carry. Use ",
         "`unit = \"record\"` for the model's own `log_lik` terms.",
         call. = FALSE)
  Y <- data$Y
  X <- data$X
  A <- data$A
  n <- nrow(Y)
  k <- ncol(Y)
  miss <- is.na(Y)
  traits <- .draw_traits(fit, NULL, k)

  eig <- if (is.null(eigen_A)) .relationship_eigen(A) else
    .check_eigen_A(eigen_A, A)
  Q <- eig$Q
  lambda <- pmax(eig$lambda, 0)
  QQ <- Q * Q

  pulled <- .individual_loo_draws(fit)
  draws <- pulled$draws
  beta <- .beta_draws(draws, k, traits, NULL, attr(fit, "x_names"))
  G <- .as_matrix_draws(.matrix_draws(draws, "G", k, traits),
                        .vector_draws(draws, "sigma_G", k, traits), k)
  E <- .as_matrix_draws(.matrix_draws(draws, "E", k, traits),
                        .vector_draws(draws, "sigma_E", k, traits), k)
  if (is.null(beta) || is.null(G) || is.null(E))
    stop("this fit does not report `beta` and the variances, so the ",
         "individual-level terms cannot be computed.", call. = FALSE)
  S <- nrow(draws)

  standard <- .standardization(Y)
  log_sd <- log(vapply(standard, `[[`, numeric(1), "sd"))
  gaps <- data.frame(individual = row(miss)[miss], trait = col(miss)[miss])
  filled <- if (nrow(gaps)) .impute_animal(draws, seq_len(S), gaps, standard)
    else NULL

  if (is.null(batch))
    batch <- max(1L, min(S, floor(2e7 / (n * k))))
  Y0 <- Y
  Y0[miss] <- 0
  out <- matrix(NA_real_, S, n)
  for (start in seq(1L, S, by = batch)) {
    take <- start:min(S, start + batch - 1L)
    R <- array(0, c(n, k, length(take)))
    for (t in seq_along(take)) {
      s <- take[t]
      value <- Y0
      if (!is.null(filled)) value[miss] <- filled[s, ]
      R[, , t] <- value - X %*% t(matrix(beta[, , s], nrow = k))
    }
    out[take, ] <- .individual_log_density(R, G[, , take, drop = FALSE],
                                           E[, , take, drop = FALSE],
                                           Q, QQ, lambda, miss)
  }
  # standardized-response density: add log sd of each observed trait
  out <- sweep(out, 2L, drop((!miss) %*% log_sd), "+")
  observed <- rowSums(!miss) > 0L
  list(log_lik = out[, observed, drop = FALSE],
       chain_id = pulled$chain_id,
       units = (rownames(Y) %||% as.character(seq_len(n)))[observed])
}

#' Leave-one-out cross-validation for a fit
#'
#' Computes the expected log pointwise predictive density with Pareto-smoothed
#' importance sampling (or WAIC), leaving out one individual at a time. Animal
#' fits with missing trait records are supported.
#'
#' @param x A fit from \code{\link{lmm_animal}} or \code{\link{lm_phen_cov}},
#'   or a draws matrix with \code{log_lik} columns.
#' @param criterion \code{"loo"} for Pareto-smoothed importance sampling, or
#'   \code{"waic"} for the widely applicable information criterion.
#' @param r_eff Relative effective sample sizes, one per observation. Computed
#'   from the within-chain correlation of the pointwise likelihoods when left
#'   \code{NULL}, which is what \code{\link[loo]{relative_eff}} recommends.
#' @param unit For animal fits, what one pointwise term is.
#'   \code{"individual"} (default) leaves out one individual: the term is the
#'   density of its observed traits given every other individual, computed after
#'   the fit from the draws of \code{G}, \code{E}, \code{beta} and, with gaps,
#'   \code{Y_mis}. \code{"record"} uses the model's own \code{log_lik}, one term
#'   per record rotated by the eigenvectors of \code{A}, which needs a complete
#'   response. Ignored for phenotypic fits and draws matrices, whose terms are
#'   already one per individual or are taken as given.
#' @param eigen_A Optional eigendecomposition of the relationship matrix the fit
#'   used (\code{attr(fit, "stanimal_data")$A}, after any reordering and
#'   \code{blend_A}), as a list with \code{Q} and \code{lambda}, to skip the
#'   \eqn{O(N^3)} decomposition for \code{unit = "individual"}. It is checked
#'   against that matrix.
#' @param ... Further arguments passed to \code{\link[loo]{loo}} or
#'   \code{\link[loo]{waic}}, such as \code{k_threshold} or \code{cores}.
#' @return An object from the \code{loo} package: a \code{psis_loo} or
#'   \code{waic} object. \code{loo_compare()} accepts several of them, which is
#'   how fits are compared; the differences in ELPD are what matters, not the
#'   absolute value. For \code{unit = "individual"}, the attribute
#'   \code{"units"} names the individual of each pointwise term, so a high
#'   Pareto \eqn{k} can be traced to an animal.
#' @details Relatives share breeding values, so an animal model does not
#'   factorize over individuals and the density of one animal depends on all the
#'   others. The default terms are those of leave-one-out for a non-factorized
#'   normal model (Bürkner, Gabry and Vehtari 2021): for each draw, the
#'   conditional density of an individual's observed traits given every other
#'   individual, with that individual's own gaps integrated out. When other
#'   individuals have gaps, the draw conditions on their sampled values
#'   (\code{Y_mis}), which keeps the importance sampling exact. Individuals with
#'   no observed trait carry no term. The cost per draw is two
#'   \eqn{N \times N} by \eqn{N \times K} matrix products, batched over draws,
#'   plus one decomposition of \code{A}.
#'
#'   Phenotypic models use one term per individual, integrating missing traits
#'   out, from their own \code{log_lik}.
#'
#'   With \code{unit = "record"}, the terms are the animal model's
#'   \code{log_lik}: for complete responses they sum to the joint log
#'   likelihood, but the rotation mixes individuals, so a term is not one animal
#'   left out. That criterion is refused for responses with gaps, whose
#'   \code{log_lik} terms use imputed values.
#'
#'   Compare fits to the same response values, observations and missingness
#'   pattern, with the same unit. The terms are densities of the standardized
#'   response, so the absolute ELPD carries a constant that depends on the
#'   data's scale; it cancels in comparisons between fits to the same response.
#'   Individual-level terms assume the built-in likelihood and are refused for a
#'   custom \code{model_file}. If response or model metadata are unavailable
#'   (a bare draws matrix or an older fit), the model's \code{log_lik} is used
#'   with a warning: the caller must verify that its terms are suitable
#'   predictive units. Extracting raw draws and passing them directly to
#'   \code{loo::loo()} or \code{loo::waic()} bypasses these checks.
#' @references Bürkner, P.-C., Gabry, J. and Vehtari, A. (2021). Efficient
#'   leave-one-out cross-validation for Bayesian non-factorized normal and
#'   Student-t models. \emph{Computational Statistics} 36, 1243--1261.
#'   \doi{10.1007/s00180-020-01045-4}
#' @seealso \code{\link{posterior_predictive}}
#' @examples
#' \dontrun{
#' fit1 <- lmm_animal(Y, X, A)
#' fit2 <- lmm_animal(Y, X, A, genetic_structure = "factor", n_factors = 2)
#' loo::loo_compare(fit_loo(fit1), fit_loo(fit2))
#' }
#' @export
fit_loo <- function(x, criterion = c("loo", "waic"), r_eff = NULL,
                    unit = c("individual", "record"), eigen_A = NULL, ...) {
  criterion <- match.arg(criterion)
  unit <- match.arg(unit)
  if (!requireNamespace("loo", quietly = TRUE))
    stop("`fit_loo()` needs the `loo` package: install.packages(\"loo\").",
         call. = FALSE)
  data <- attr(x, "stanimal_data")
  individual <- unit == "individual" && is.list(data) &&
    identical(data$family, "animal")
  if (individual) {
    # the refusals for Laplace and Gibbs objects live in .log_lik_draws()
    if (inherits(x, "stanimal_laplace") || .is_gibbs_draws(x))
      .log_lik_draws(x)
    draws <- .individual_log_lik(x, eigen_A = eigen_A)
  } else {
    draws <- .log_lik_draws(x)
  }
  log_lik <- draws$log_lik
  result <- if (criterion == "waic") {
    loo::waic(log_lik, ...)
  } else {
    if (is.null(r_eff))
      r_eff <- loo::relative_eff(exp(log_lik), chain_id = draws$chain_id)
    loo::loo(log_lik, r_eff = r_eff, ...)
  }
  if (individual) attr(result, "units") <- draws$units
  result
}

#' @rdname fit_loo
#' @param x A fit object.
#' @method loo CmdStanMCMC
#' @keywords internal
loo.CmdStanMCMC <- function(x, ...) fit_loo(x, ...)

#' @rdname fit_loo
#' @method loo stanfit
#' @keywords internal
loo.stanfit <- function(x, ...) fit_loo(x, ...)

# `loo` is only suggested, so its methods cannot go through `S3method()` in
# NAMESPACE: R would look for the generic at load time and fail. Register them
# instead, now if `loo` is already loaded and otherwise when it loads, which is
# what makes `loo(fit)` work without an import.
.onLoad <- function(libname, pkgname) {
  register_loo_methods <- function() {
    registerS3method("loo", "CmdStanMCMC", loo.CmdStanMCMC,
                     envir = asNamespace("loo"))
    registerS3method("loo", "stanfit", loo.stanfit, envir = asNamespace("loo"))
  }
  setHook(packageEvent("loo", "onLoad"), function(...) register_loo_methods())
  if (isNamespaceLoaded("loo"))
    register_loo_methods()
}
