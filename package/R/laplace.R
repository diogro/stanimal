#' Laplace approximations of the covariance fits
#'
#' Fits the same model as \code{\link{lmm_animal}} or \code{\link{lm_phen_cov}}
#' and returns a Gaussian approximation at the posterior mode instead of
#' posterior draws. It can be much faster than sampling, but has no MCMC
#' diagnostics: there is no \code{Rhat}, divergence count or effective sample
#' size, because there are no chains.
#'
#' @param draws Number of draws from the approximation used for its reported
#'   means, standard deviations and median absolute deviations. A whole number.
#' @param seed Optional seed for optimization and approximation draws. If omitted,
#'   a seed is drawn from R's random-number generator and recorded in the result.
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
#'       original data scale. Quantities are labelled with the column names of
#'       \code{Y} (and \code{X} for \code{beta}), e.g. \code{G[wing,tail]}
#'       and \code{h2[tail]}, as \code{fit_draws()} names sampler fits; an
#'       unnamed \code{Y} keeps Stan's indices (\code{G[2,1]}).
#'     \item \code{draws} the draws from the approximation, for inspection
#'       only, with the same labels.
#'     \item \code{y_names}, \code{x_names} the column names of \code{Y} and
#'       \code{X}.
#'     \item \code{elapsed} seconds the fit took.
#'     \item \code{call} the matched call.
#'   }
#'   It is not a posterior sample and the draw-based functions refuse it. Use it
#'   to check that the sampler and a Laplace approximation find the same region,
#'   to screen many trait sets quickly, or to get a covariance when a full fit is
#'   too slow, and not to report a heritability with an interval.
#'
#' @details The approximation is Gaussian on the unconstrained parameter scale,
#'   centred on the mode, with covariance taken from the curvature there.
#'   Quantities on their reported scales are transformations of those draws.
#'   The approximation is poor wherever the
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
#'   Optimization and CmdStan failures and warnings are surfaced; there is no
#'   automatic retry. If the Laplace stage fails, its CmdStan output is included
#'   in the error. The mode is written with high precision to avoid loss of
#'   precision when CmdStan reads a correlation Cholesky factor back in.
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
  .warn_dots(...)
  beta_prior_supplied <- !missing(beta_prior)
  lkj_prior_supplied <- !missing(lkj_prior)
  parametrization_supplied <- !missing(parametrization)
  hs_supplied <- c(hs_global_scale = !missing(hs_global_scale),
                   hs_slab_scale = !missing(hs_slab_scale),
                   column_decay = !missing(column_decay))
  parametrization <- match.arg(parametrization)
  residual_structure <- match.arg(residual_structure)
  genetic_structure <- match.arg(genetic_structure)
  for (arg in names(hs_supplied))
    .warn_ignored(hs_supplied[[arg]], genetic_structure == "factor", arg)
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
  .warn_ignored(beta_prior_supplied,
                (ncol(Y) == 1L && parametrization == "h2") || factor_model ||
                  h2_param, "beta_prior")
  .warn_ignored(lkj_prior_supplied,
                .animal_lkj_used(ncol(Y), genetic_structure, residual_structure),
                "lkj_prior")
  .warn_ignored(parametrization_supplied, !factor_model, "parametrization")

  relationship <- .animal_relationship(A, Y, blend_A = blend_A,
                                       need = "for the Laplace approximation")
  A <- relationship$A

  if (ncol(Y) == 1L) {
    stan_data <- list(J = ncol(X), N = nrow(Y), X = X, Y = as.numeric(Y),
                      Q = relationship$Q,
                      lambda = relationship$lambda, miss = as.integer(missing$miss),
                      n_missing = missing$n_missing,
                      mis_row = as.array(missing$mis_row),
                      x_mean_weight = YX$x_mean_weight)
    # a custom model may still declare A; extra data are ignored
    if (!is.null(model_file)) stan_data$A <- A
    model_code <- .read_animal_model_code(
      model_file,
      if (parametrization == "h2") "animalModelUniH2.stan" else "animalModelUni.stan"
    )
    if (parametrization == "h2") stan_data$beta_prior <- beta_prior
    quantities <- c("sigma_G", "sigma_E", "h2")
  } else {
    stan_data <- list(K = ncol(Y), J = ncol(X), N = nrow(Y), X = X, Y = Y,
                      Q = relationship$Q,
                      lambda = relationship$lambda, miss = missing$miss,
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
    # a custom model may still declare A; extra data are ignored
    if (!is.null(model_file)) stan_data$A <- A
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

  .laplace_fit(model_code, stan_data, quantities, draws, seed, call,
               y_names = colnames(Y), x_names = colnames(X))
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
  .warn_dots(...)
  lkj_prior_supplied <- !missing(lkj_prior)
  hs_supplied <- c(hs_global_scale = !missing(hs_global_scale),
                   hs_slab_scale = !missing(hs_slab_scale),
                   column_decay = !missing(column_decay))
  covariance_structure <- match.arg(covariance_structure)
  for (arg in names(hs_supplied))
    .warn_ignored(hs_supplied[[arg]], covariance_structure == "factor", arg)
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
  .warn_ignored(lkj_prior_supplied, !factor_model, "lkj_prior")

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

  .laplace_fit(model_code, stan_data, quantities, draws, seed, call,
               y_names = colnames(Y), x_names = colnames(X))
}

# `...` is kept for symmetry with the samplers but nothing consumes it here
.warn_dots <- function(...) {
  if (...length() > 0L) {
    arg_names <- names(list(...))
    if (is.null(arg_names)) arg_names <- rep("", ...length())
    arg_names[!nzchar(arg_names)] <- "<unnamed>"
    warning("Arguments in `...` are ignored by the Laplace approximation: ",
            paste(arg_names, collapse = ", "), ".",
            call. = FALSE)
  }
  invisible(NULL)
}

# Run the Laplace approximation on a compiled model and wrap the summaries.
.laplace_fit <- function(model_code, stan_data, quantities, draws, seed, call,
                         y_names = NULL, x_names = NULL) {
  draws <- .check_count(draws, "draws")
  if (!is.null(seed)) seed <- .check_whole(seed, "seed")
  if (!requireNamespace("cmdstanr", quietly = TRUE))
    stop("A Laplace approximation needs the cmdstanr package; install it, or ",
         "use the sampler with `lmm_animal()` or `lm_phen_cov()`.", call. = FALSE)

  stan_file <- cmdstanr::write_stan_file(model_code, dir = stan_cache_dir())
  mod <- cmdstanr::cmdstan_model(stan_file)

  seed <- if (is.null(seed)) sample.int(1e6, 1L) else seed
  start <- proc.time()[["elapsed"]]
  result <- .laplace_run(mod, stan_data, quantities, draws, seed)
  elapsed <- proc.time()[["elapsed"]] - start

  summary <- result$summary
  colnames(summary) <- c("quantity", "mean", "sd", "mad")
  # Label trait and predictor indices with the names of Y and X, as
  # fit_draws() does for sampler fits: G[2,1] -> G[wing,tail].
  summary$quantity <- .label_quantities(summary$quantity, y_names, x_names)
  laplace_draws <- result$fit$draws()
  posterior::variables(laplace_draws) <- .label_quantities(
    posterior::variables(laplace_draws), y_names, x_names)

  structure(
    list(summary = summary,
         draws = laplace_draws,
         y_names = y_names,
         x_names = x_names,
         elapsed = elapsed,
         draws_requested = draws,
         seed = result$seed,
         call = call),
    class = "stanimal_laplace"
  )
}

.laplace_run <- function(model, stan_data, quantities, draws, seed) {
  # No unique output_basename here (see .unique_output()): released cmdstanr
  # (0.9.0) reuses it for the internal mode file, and CmdStan then refuses
  # the run because the mode and output files share a name. Nor is one needed:
  # the summary and draws are read below, before this function returns, and
  # the fit itself is not kept, so a later run cannot change them.
  fit <- model$laplace(data = stan_data, draws = draws, refresh = 0,
                       seed = seed, sig_figs = 18, show_messages = FALSE)
  codes <- fit$return_codes()
  if (!length(codes) || anyNA(codes) || any(codes != 0)) {
    output <- utils::capture.output(fit$output())
    stop("The Laplace approximation failed (seed ", seed, "). CmdStan output:\n",
         paste(output, collapse = "\n"), call. = FALSE)
  }
  summary <- fit$summary(quantities)[, c("variable", "mean", "sd", "mad")]
  list(fit = fit, summary = summary, seed = seed)
}

#' @export
print.stanimal_laplace <- function(x, digits = 3, ...) {
  cat("Laplace approximation, not a posterior sample\n")
  cat("  model :", deparse(x$call[[1L]]), "\n")
  cat("  draws :", x$draws_requested, " (seed ", x$seed, ")\n", sep = "")
  cat("  time  :", round(x$elapsed, 1), "s\n")
  cat("\nNo MCMC convergence diagnostics exist for this object, and the draw-based\n",
      "functions refuse it. Standard deviations for variances and\n",
      "heritabilities are the least reliable entries.\n\n", sep = "")
  print(as.data.frame(x$summary[, c("quantity", "mean", "sd")]),
        digits = digits, row.names = FALSE)
  invisible(x)
}

# Replace numeric trait (and predictor) indices in Stan variable names with
# the names of Y (and X). Trait-by-trait matrices get both indices labelled,
# per-trait vectors their one index, loadings their trait index (the factor
# index stays a number), and beta its trait and predictor. Anything else, and
# everything when there are no names, is returned unchanged.
.label_quantities <- function(vars, y_names = NULL, x_names = NULL) {
  if (is.null(y_names)) return(vars)
  by_trait_2 <- c("G", "E", "P", "corrG", "corrE", "corrP")
  by_trait_1 <- c("h2", "sigma_G", "sigma_E", "sigma_P", "s_G", "uniqueness",
                  "communality", "genetic_uniqueness", "residual_uniqueness")
  m <- regmatches(vars, regexec("^([A-Za-z_][A-Za-z0-9_]*)\\[([0-9]+)(,([0-9]+))?\\]$",
                                vars))
  out <- vars
  for (k in seq_along(vars)) {
    p <- m[[k]]
    if (!length(p)) next
    base <- p[2L]
    i <- as.integer(p[3L])
    j <- if (nzchar(p[5L])) as.integer(p[5L]) else NA_integer_
    label <- function(names, n) if (!is.na(n) && n <= length(names)) names[n] else NA
    new <- if (base %in% by_trait_2 && !is.na(j)) {
      c(label(y_names, i), label(y_names, j))
    } else if (base %in% by_trait_1 && is.na(j)) {
      label(y_names, i)
    } else if (base == "loadings" && !is.na(j)) {
      c(label(y_names, i), j)
    } else if (base == "beta" && !is.na(j) && !is.null(x_names)) {
      c(label(y_names, i), label(x_names, j))
    } else NULL
    if (!is.null(new) && !anyNA(new))
      out[k] <- paste0(base, "[", paste(new, collapse = ","), "]")
  }
  out
}
