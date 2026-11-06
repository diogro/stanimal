# Predictive model comparison for a fitted model.
#
# Suitable fits report one `log_lik` term per record. Keep the chain structure
# intact when handing those draws to `loo` for its effective sample size correction.
#
# The terms are on the rotated scale the model block works on, which is where
# the records are independent for complete animal responses. One term is one
# *rotated* record, not one individual. Missing animal responses are imputed
# parameters, so those terms are not suitable for an observed-data criterion.

.check_loo_data <- function(fit) {
  data <- attr(fit, "stanimal_data")
  if (is.null(data) || !is.list(data) ||
      !(identical(data$family, "animal") ||
        identical(data$family, "phenotypic")) || is.null(data$Y)) {
    warning("Cannot verify whether `log_lik` is suitable for LOO or WAIC: ",
            "response and model metadata are unavailable. Ensure the terms ",
            "are appropriate observed-data predictive units; incomplete animal ",
            "fits' completed-data terms are not.", call. = FALSE)
  } else if (identical(data$family, "animal") && anyNA(data$Y)) {
    stop("LOO and WAIC are unavailable for animal fits with missing responses: ",
         "their `log_lik` uses imputed values rather than an observed-data ",
         "likelihood. Use `posterior_predictive()` for a model check.",
         call. = FALSE)
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

#' Leave-one-out cross-validation for a fit
#'
#' Computes the expected log pointwise predictive density with Pareto-smoothed
#' importance sampling, by handing the model's `log_lik` draws to
#' \code{\link[loo]{loo}}. Built-in models report that quantity, but incomplete
#' animal fits cannot be compared this way (see Details).
#'
#' @param x A fit from \code{\link{lmm_animal}} or \code{\link{lm_phen_cov}},
#'   or a draws matrix with \code{log_lik} columns.
#' @param criterion \code{"loo"} for Pareto-smoothed importance sampling, or
#'   \code{"waic"} for the widely applicable information criterion.
#' @param r_eff Relative effective sample sizes, one per observation. Computed
#'   from the within-chain correlation of the pointwise likelihoods when left
#'   \code{NULL}, which is what \code{\link[loo]{relative_eff}} recommends.
#' @param ... Further arguments passed to \code{\link[loo]{loo}} or
#'   \code{\link[loo]{waic}}, such as \code{k_threshold} or \code{cores}.
#' @return An object from the \code{loo} package: a \code{psis_loo} or
#'   \code{waic} object. \code{loo_compare()} accepts several of them, which is
#'   how fits are compared; the differences in ELPD are what matters, not the
#'   absolute value.
#' @details Phenotypic models use one term per individual, integrating missing
#'   traits out. Animal models use one term per rotated response record. For
#'   complete animal responses the terms sum to the joint log likelihood, but
#'   the rotation can mix individuals, so a term is not one animal left out.
#'   Read the result as a comparison between
#'   model structures or fixed effects on the same response, observations and
#'   missingness pattern, with corresponding predictive units. Animal fits
#'   must also have the same relationship matrix and eigenbasis. Different
#'   response traits or values are not directly comparable. This is not an
#'   individual-level predictive check.
#'
#'   The terms are the density of the standardized response, so the absolute
#'   ELPD carries a constant that depends on the data's scale; it cancels in
#'   comparisons between fits to the same response.
#'   Animal fits with missing responses are refused, including custom-model
#'   fits: the built-in terms use imputed values, not the marginalized
#'   observed-data likelihood. \code{\link{posterior_predictive}} remains
#'   available for model checking. If response or model metadata are unavailable
#'   (including for a bare draws matrix or older fit), a warning is issued:
#'   the caller must verify that the terms are suitable predictive units.
#'   Custom models may define different likelihoods, so their suitability
#'   cannot be verified from metadata alone. Extracting raw draws and passing
#'   them directly to \code{loo::loo()} or \code{loo::waic()} bypasses these
#'   checks.
#' @seealso \code{\link{posterior_predictive}}
#' @examples
#' \dontrun{
#' fit1 <- lmm_animal(Y, X, A)
#' fit2 <- lmm_animal(Y, X, A, genetic_structure = "factor", n_factors = 2)
#' loo::loo_compare(fit_loo(fit1), fit_loo(fit2))
#' }
#' @export
fit_loo <- function(x, criterion = c("loo", "waic"), r_eff = NULL, ...) {
  criterion <- match.arg(criterion)
  if (!requireNamespace("loo", quietly = TRUE))
    stop("`fit_loo()` needs the `loo` package: install.packages(\"loo\").",
         call. = FALSE)
  draws <- .log_lik_draws(x)
  log_lik <- draws$log_lik
  if (criterion == "waic")
    return(loo::waic(log_lik, ...))
  if (is.null(r_eff))
    r_eff <- loo::relative_eff(exp(log_lik), chain_id = draws$chain_id)
  loo::loo(log_lik, r_eff = r_eff, ...)
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
