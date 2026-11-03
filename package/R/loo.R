# Predictive model comparison for a fitted model.
#
# Every model in the package reports one pointwise `log_lik` term per record, so
# leave-one-out cross-validation is a matter of handing those draws to the `loo`
# package with the chain structure intact (the effective sample size correction
# needs to know which draws belong to which chain).
#
# The terms are on the rotated scale the model block works on, which is where
# the records are independent. That is what makes the sum exact, and it is also
# the caveat: one term is one *rotated* record, not one individual, so this
# compares model structures on the same data rather than leaving out an animal.

# Draws of `log_lik` with their chain structure: a draws by records matrix, as
# `loo` wants, and the chain each draw came from.
.log_lik_draws <- function(fit) {
  if (.is_gibbs_draws(fit))
    stop("`lmm_animal_gibbs()` does not report `log_lik`, so LOO and WAIC ",
         "cannot be computed for it. Use `lmm_animal()` with `cmd_stan = FALSE` ",
         "for a Stan fit of the same model.", call. = FALSE)
  if (inherits(fit, "stanfit")) {
    draws <- tryCatch(as.array(fit, pars = "log_lik"),
                      error = function(e) NULL)
  } else if (inherits(fit, "CmdStanFit") ||
             (is.list(fit) && is.function(fit$draws))) {
    draws <- tryCatch(
      as.array(fit$draws(variables = "log_lik", format = "draws_array")),
      error = function(e) NULL)
  } else if (is.matrix(fit) && !is.null(colnames(fit))) {
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
#' \code{\link[loo]{loo}}. Every model in the package reports that quantity, so
#' this works for univariate, multivariate, factor and phenotypic fits alike.
#'
#' @param x A fit from \code{\link{lmm_animal}} or \code{\link{lm_phen_cov}}.
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
#' @details The pointwise terms are those of the model block, one per record of
#'   the rotated response the likelihood factorizes over. When \code{A} is the
#'   identity that is one term per individual; otherwise the rotation mixes
#'   individuals, so the terms still sum to the exact joint log likelihood but a
#'   term is not one animal left out. Read the result as a comparison between
#'   model structures on the same data — unstructured against factor, different
#'   trait sets, different fixed effects — and not as an individual-level
#'   predictive check. For the latter, use \code{\link{posterior_predictive}}.
#'
#'   The terms are the density of the standardized response, so the absolute
#'   ELPD carries a constant that depends on the data's scale; it cancels in
#'   comparisons between fits to the same response.
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
