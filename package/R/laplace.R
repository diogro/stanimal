#' Laplace approximations of the covariance fits
#'
#' Fits the same model as \code{\link{lmm_animal}} or \code{\link{lm_phen_cov}}
#' and returns a Gaussian approximation at the posterior mode instead of
#' posterior draws. It is much faster than sampling and it comes with no
#' diagnostics: there is no \code{Rhat}, no divergence count and no effective
#' sample size, because there are no chains.
#'
#' @param draws Number of draws to take from the approximation, used only for the
#'   posterior standard deviations. A whole number.
#' @param seed Optional seed for the draws.
#' @inheritParams lmm_animal
#' @param Y,X,A Relationship matrix and data, as in \code{\link{lmm_animal}}.
#'   \code{lmm_animal_laplace()} needs \code{A}; \code{lm_phen_cov_laplace()}
#'   does not take one.
#' @param covariance_structure Covariance structure for
#'   \code{\link{lm_phen_cov_laplace}}, as in \code{\link{lm_phen_cov}}:
#'   \code{"unstructured"} (default) or \code{"factor"}.
#' @return An object of class \code{stanimal_laplace}, a list with
#'   \itemize{
#'     \item \code{summary} a data frame of posterior means and standard
#'       deviations for the covariance quantities the model reports, on the
#'       original data scale.
#'     \item \code{draws} the draws from the approximation, for inspection only.
#'     \item \code{elapsed} seconds the fit took.
#'     \item \code{call} the matched call.
#'   }
#'   It is not a posterior sample and the draw-based functions refuse it. Use it
#'   to check that the sampler and a Laplace approximation find the same region,
#'   to screen many trait sets quickly, or to get a covariance when a full fit is
#'   too slow, and not to report a heritability with an interval.
#'
#' @details The approximation is a Gaussian centred on the posterior mode, with
#'   the covariance taken from the curvature there. It is poor wherever the
#'   posterior is skewed, and the diagonal of a covariance matrix is the first
#'   place that happens, so a standard deviation for a variance or a heritability
#'   is the least trustworthy number these functions return. In bivariate animal
#'   models with 120 individuals the means came out systematically below the
#'   sampler's, by up to 0.18 on entries of order one, because the mode sits
#'   below the mean when the posterior is skewed, and the standard deviation of
#'   the smallest variance was about 40\% too wide.
#'
#'   Only \code{cmdstanr} can run the approximation; there is no \code{rstan}
#'   equivalent, so the functions error when it is not installed rather than
#'   falling back to sampling.
#'
#' @seealso \code{\link{lmm_animal}}, \code{\link{lm_phen_cov}}
#' @name laplace
NULL

#' @rdname laplace
#' @export
lmm_animal_laplace <- function(Y, X = NULL, A,
                               parametrization = c("variance", "h2"),
                               lkj_prior = 2, beta_prior = c(1, 1),
                               blend_A = FALSE, model_file = NULL,
                               residual_structure = c("unstructured",
                                                      "diagonal",
                                                      "shared_factor"),
                               genetic_structure = c("unstructured", "factor"),
                               n_factors = NULL, hs_global_scale = 1,
                               hs_slab_scale = 1, column_decay = 0.8,
                               draws = 2000, seed = NULL, ...) {
  call <- match.call()
  parametrization <- match.arg(parametrization)
  residual_structure <- match.arg(residual_structure)
  genetic_structure <- match.arg(genetic_structure)
  lkj_prior <- .check_positive(lkj_prior, "lkj_prior")
  beta_prior <- .check_beta_prior(beta_prior)
  hs_global_scale <- .check_positive(hs_global_scale, "hs_global_scale")
  hs_slab_scale <- .check_positive(hs_slab_scale, "hs_slab_scale")
  column_decay <- .check_in_range(column_decay, "column_decay", 0, 1,
                                  open_low = TRUE)
  blend_A <- .check_flag(blend_A, "blend_A")

  YX <- .check_response_design(Y, X, allow_na = TRUE, allow_empty_rows = TRUE)
  Y <- YX$Y
  X <- YX$X
  missing <- YX$missing
  A <- as.matrix(A)

  factor_model <- genetic_structure == "factor"
  if (residual_structure == "shared_factor" && !factor_model)
    stop("`residual_structure = \"shared_factor\"` requires ",
         "`genetic_structure = \"factor\"`.", call. = FALSE)
  shared_factor_model <- factor_model &&
    residual_structure == "shared_factor"
  h2_param <- ncol(Y) > 1L && !factor_model && parametrization == "h2"

  if (factor_model) {
    if (ncol(Y) < 2L)
      stop("`genetic_structure = \"factor\"` needs at least two traits.",
           call. = FALSE)
    if (!is.null(model_file))
      stop("`model_file` cannot be combined with ",
           "`genetic_structure = \"factor\"`.", call. = FALSE)
    if (is.null(n_factors)) {
      n_factors <- .default_n_factors(ncol(Y))
      warning("`n_factors` was not supplied; using ", n_factors,
              " factors for ", ncol(Y), " traits.", call. = FALSE)
    }
    n_factors <- .check_count(n_factors, "n_factors")
    if (n_factors >= ncol(Y))
      stop("`n_factors` must be smaller than the number of traits (",
           ncol(Y), ").", call. = FALSE)
  } else if (!is.null(n_factors)) {
    stop("`n_factors` is only used with `genetic_structure = \"factor\"`.",
         call. = FALSE)
  }

  n <- nrow(Y)
  A <- .check_relationship_matrix(A, Y)
  A <- .check_symmetric(A)
  eigenvalues <- .check_positive_semidefinite(A)
  if (blend_A) {
    A <- 0.99 * A + 0.01 * diag(n)
    eigenvalues <- .relationship_eigenvalues(A)
  }
  if (min(eigenvalues) <= .eigen_tolerance(eigenvalues))
    stop("`A` must be positive definite for the Laplace approximation; set ",
         "`blend_A = TRUE` to regularize it.", call. = FALSE)

  if (ncol(Y) == 1L) {
    stan_data <- list(J = ncol(X), N = nrow(Y), X = X, Y = as.numeric(Y),
                      A = A, miss = as.integer(missing$miss),
                      n_missing = missing$n_missing,
                      mis_row = as.array(missing$mis_row),
                      x_mean_weight = YX$x_mean_weight)
    model_code <- .read_animal_model_code(
      model_file,
      if (parametrization == "h2") "animalModelUniH2.stan" else "animalModelUni.stan"
    )
    if (parametrization == "h2") stan_data$beta_prior <- beta_prior
    quantities <- c("sigma_G", "sigma_E", "h2")
  } else {
    stan_data <- list(K = ncol(Y), J = ncol(X), N = nrow(Y), X = X, Y = Y,
                      A = A, miss = missing$miss,
                      n_missing = missing$n_missing,
                      mis_row = as.array(missing$mis_row),
                      mis_col = as.array(missing$mis_col),
                      x_mean_weight = YX$x_mean_weight)
    if (shared_factor_model) {
      stan_data$M <- n_factors
      stan_data$hs_global_scale <- hs_global_scale
      stan_data$hs_slab_scale <- hs_slab_scale
      stan_data$column_decay <- column_decay
      stan_data$h2_prior_a <- beta_prior[1L]
      stan_data$h2_prior_b <- beta_prior[2L]
    } else {
      stan_data$lkj_prior <- lkj_prior
      stan_data$unstructured_E <-
        as.integer(residual_structure == "unstructured")
      if (h2_param) stan_data$beta_prior <- beta_prior
      if (factor_model) {
        stan_data$M <- n_factors
        stan_data$hs_global_scale <- hs_global_scale
        stan_data$hs_slab_scale <- hs_slab_scale
        stan_data$column_decay <- column_decay
        stan_data$beta_prior <- beta_prior
      }
    }
    default_model <- if (shared_factor_model) {
      "animalModelGSharedFactor.stan"
    } else if (factor_model) {
      "animalModelGFactor.stan"
    } else if (h2_param) {
      "animalModelGH2.stan"
    } else {
      "animalModelG.stan"
    }
    model_code <- .read_animal_model_code(model_file, default_model)
    quantities <- c("G", "E", "P", "corrG", "sigma_G", "sigma_E", "h2")
    if (shared_factor_model) {
      quantities <- c(quantities, "loadings", "uniqueness", "communality",
                      "factor_magnitude", "factor_heritability",
                      "genetic_uniqueness", "residual_uniqueness")
    } else if (factor_model) {
      quantities <- c(quantities, "loadings", "uniqueness", "communality",
                      "factor_magnitude")
    }
  }

  .laplace_fit(model_code, stan_data, quantities, draws, seed, call)
}

#' @rdname laplace
#' @export
lm_phen_cov_laplace <- function(Y, X = NULL, lkj_prior = 2,
                                covariance_structure = c("unstructured",
                                                         "factor"),
                                n_factors = NULL, hs_global_scale = 1,
                                hs_slab_scale = 1, column_decay = 0.8,
                                model_file = NULL, draws = 2000, seed = NULL,
                                ...) {
  call <- match.call()
  covariance_structure <- match.arg(covariance_structure)
  lkj_prior <- .check_positive(lkj_prior, "lkj_prior")
  hs_global_scale <- .check_positive(hs_global_scale, "hs_global_scale")
  hs_slab_scale <- .check_positive(hs_slab_scale, "hs_slab_scale")
  column_decay <- .check_in_range(column_decay, "column_decay", 0, 1,
                                  open_low = TRUE)

  YX <- .check_response_design(Y, X, allow_na = TRUE)
  Y <- YX$Y
  X <- YX$X
  missing <- YX$missing
  factor_model <- covariance_structure == "factor"

  if (factor_model) {
    if (ncol(Y) < 2L)
      stop("`covariance_structure = \"factor\"` needs at least two traits.",
           call. = FALSE)
    if (!is.null(model_file))
      stop("`model_file` cannot be combined with ",
           "`covariance_structure = \"factor\"`.", call. = FALSE)
    if (is.null(n_factors)) {
      n_factors <- .default_n_factors(ncol(Y))
      warning("`n_factors` was not supplied; using ", n_factors,
              " factors for ", ncol(Y), " traits.", call. = FALSE)
    }
    n_factors <- .check_count(n_factors, "n_factors")
    if (n_factors >= ncol(Y))
      stop("`n_factors` must be smaller than the number of traits (",
           ncol(Y), ").", call. = FALSE)
  } else if (!is.null(n_factors)) {
    stop("`n_factors` is only used with `covariance_structure = \"factor\"`.",
         call. = FALSE)
  }

  stan_data <- list(K = ncol(Y), J = ncol(X), N = nrow(Y), X = X, Y = Y,
                    miss = missing$miss, n_pattern = missing$n_pattern,
                    pattern_obs = missing$pattern_obs,
                    row_order = as.array(missing$row_order),
                    n_rows = as.array(missing$n_rows),
                    first_row = as.array(missing$first_row),
                    x_mean_weight = YX$x_mean_weight)
  if (factor_model) {
    stan_data$M <- n_factors
    stan_data$hs_global_scale <- hs_global_scale
    stan_data$hs_slab_scale <- hs_slab_scale
    stan_data$column_decay <- column_decay
  } else {
    stan_data$lkj_prior <- lkj_prior
  }
  model_code <- .read_animal_model_code(
    model_file,
    if (factor_model) "phenotypicModelFactor.stan" else "phenotypicModel.stan"
  )
  quantities <- c("P", "corrP", "sigma_P")
  if (factor_model)
    quantities <- c(quantities, "loadings", "uniqueness", "communality",
                    "factor_magnitude")

  .laplace_fit(model_code, stan_data, quantities, draws, seed, call)
}

# Run the Laplace approximation on a compiled model and wrap the summaries.
.laplace_fit <- function(model_code, stan_data, quantities, draws, seed, call,
                         attempts = 3L) {
  if (!requireNamespace("cmdstanr", quietly = TRUE))
    stop("A Laplace approximation needs the cmdstanr package; install it, or ",
         "use the sampler with `lmm_animal()` or `lm_phen_cov()`.", call. = FALSE)
  draws <- .check_count(draws, "draws")
  if (!is.null(seed)) seed <- .check_whole(seed, "seed")

  stan_file <- cmdstanr::write_stan_file(model_code, dir = stan_cache_dir())
  mod <- cmdstanr::cmdstan_model(stan_file)

  # A draw from the approximation sometimes fails the exact unit-vector check
  # when CmdStan unconstrains a Cholesky factor of a correlation matrix, which
  # aborts that fit. It is a property of the seed rather than of the model, so
  # the seed is advanced and the fit retried before giving up.
  seeds <- if (is.null(seed)) sample.int(1e6, attempts) else
    seed + seq_len(attempts) - 1L
  result <- NULL
  elapsed <- 0
  for (s in seeds) {
    args <- list(data = stan_data, draws = draws, refresh = 0, seed = s)
    start <- proc.time()[["elapsed"]]
    # the optimiser prints its convergence line to stdout, which is noise inside
    # a function that is meant to be called from a script
    value <- tryCatch({
      invisible(utils::capture.output(fit <- do.call(mod$laplace, args)))
      summary <- fit$summary(quantities)[, c("variable", "mean", "sd", "mad")]
      list(fit = fit, summary = summary, seed = s)
    }, error = function(e) NULL)
    elapsed <- elapsed + proc.time()[["elapsed"]] - start
    if (!is.null(value)) {
      result <- value
      break
    }
  }
  if (is.null(result))
    stop("The Laplace approximation failed on ", attempts, " seeds. CmdStan ",
         "rejects a draw when a Cholesky factor of a correlation matrix is not ",
         "a unit vector to machine precision, and that check can fail for a ",
         "given mode. Try `draws` larger, a different `seed`, or sample the ",
         "model with `lmm_animal()` / `lm_phen_cov()`.", call. = FALSE)

  summary <- result$summary
  colnames(summary) <- c("quantity", "mean", "sd", "mad")

  structure(
    list(summary = summary,
         draws = result$fit$draws(),
         elapsed = elapsed,
         draws_requested = draws,
         seed = result$seed,
         call = call),
    class = "stanimal_laplace"
  )
}

#' @export
print.stanimal_laplace <- function(x, digits = 3, ...) {
  cat("Laplace approximation, not a posterior sample\n")
  cat("  model :", deparse(x$call[[1L]]), "\n")
  cat("  draws :", x$draws_requested, " (seed ", x$seed, ")\n", sep = "")
  cat("  time  :", round(x$elapsed, 1), "s\n")
  cat("\nNo convergence diagnostics exist for this object, and the draw-based\n",
      "functions refuse it. Standard deviations for variances and\n",
      "heritabilities are the least reliable entries.\n\n", sep = "")
  print(as.data.frame(x$summary[, c("quantity", "mean", "sd")]),
        digits = digits, row.names = FALSE)
  invisible(x)
}
