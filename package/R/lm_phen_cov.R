#' Fit a Bayesian phenotypic covariance matrix
#'
#' Estimates the phenotypic covariance matrix \code{P} of one or more traits
#' after accounting for fixed effects, for data with no relatedness or genomic
#' information. It is the counterpart of \code{\link{lmm_animal}} without the
#' genetic component, and uses the same priors, scaling and outputs.
#' @param Y Numeric response vector or matrix, with individuals in rows. Every trait must vary: constant traits cannot be standardized and are an error. Missing values (\code{NA}) are allowed and are not discarded: see Details.
#' @param X Optional finite numeric fixed-effects matrix with the same rows as Y and full column rank. The default \code{NULL} fits an intercept-only model: one column of ones, named \code{"(Intercept)"}. A matrix must be able to represent a constant (an intercept column, or the 0/1 indicators of every level of a factor), so that the mean of the response can be reported. Continuous columns are centered and scaled internally; returned beta draws use original X units.
#' @param lkj_prior Shape parameter of the LKJ prior on the correlation matrix. \code{1} is uniform; values above \code{1} shrink correlations toward zero and values below \code{1} favor larger correlations. The default \code{2} is recommended for many traits, because at high dimensionality the uniform prior concentrates mass near \code{+/-1} and overestimates the magnitude of correlations. A finite number above 0. It is used only by \code{covariance_structure = "unstructured"}; explicitly supplying it for the factor model warns that it is ignored.
#' @param covariance_structure Covariance structure for the phenotypic model: "unstructured" (default, scales times an LKJ correlation matrix) or "factor", which uses \eqn{P = \mathrm{diag}(\sigma_P)(\Lambda\Lambda^T + \Psi)\mathrm{diag}(\sigma_P)} with regularized-horseshoe loadings, where \eqn{\sigma_P} are the residual standard deviations, estimated with the same half-normal prior as the unstructured model.
#' @param n_factors Number of factors when \code{covariance_structure = "factor"}, a whole number between 1 and the number of traits minus one. If \code{NULL}, the model falls back to the same sublinear rule as \code{lmm_animal()} and warns with the value; that fallback depends only on the number of traits, so use \code{\link{suggest_n_factors}} to choose it from the data when the factors themselves are meant to be interpreted. Ignored otherwise.
#' @param hs_global_scale,hs_slab_scale Scales of the regularized-horseshoe prior on the factor loadings. Both must be positive and are only used by the factor model.
#' @param column_decay Standard-deviation multiplier for successive factor columns, between 0 and 1. The default 0.8 orders the columns so that the factor labels stay meaningful; with \code{column_decay = 1} the columns are exchangeable a priori and their labels can swap between chains. Few factors need a stronger decay than many: with two columns of comparable size 0.8 is often too mild to separate them, while \code{0.5^(m - 1)} would shrink the fifth column of a many-factor model to a sixteenth of the first. Only used by the factor model, and explicitly supplying it otherwise warns that it is ignored.
#' @param iter Number of MCMC iterations, a whole number.
#' @param warmup Number of warmup iterations, a whole number smaller than \code{iter}.
#' @param chains Number of chains to run, a whole number.
#' @param cores Number of chains to run in parallel (cmdstanr only), a whole number.
#' @param control Named list of sampler controls, passed to whichever engine runs: \code{adapt_delta}, \code{adapt_engaged}, \code{max_treedepth}, \code{step_size} (rstan's \code{stepsize} is accepted), \code{init_buffer}, \code{term_buffer}, \code{window} and \code{metric}. Any other name is an error, so a control is never silently ignored.
#' @param cmd_stan Logical, use cmdstanr. If FALSE, uses rstan. Both engines are suggested packages: when cmdstanr is not installed the fit falls back to rstan with a warning, and when neither is installed the fit errors.
#' @param model_file Optional path to a custom Stan file replacing the default model.
#' @param ... additional parameters for stan sampling.
#' @return With \code{cmd_stan = TRUE}, a \code{CmdStanMCMC} object from cmdstanr, usable with posterior and bayesplot (see \code{\link{lmm_animal}}); use \code{\link{rename_beta}} for \code{beta} draws named after the columns of \code{X}, and \code{\link{fit_draws}} or \code{\link{draws_mean}} to get \code{P}, \code{corrP} and \code{sigma_P} as matrices. With \code{cmd_stan = FALSE}, a \code{stanfit} object.
#'
#' The unstructured fit reports \code{P} (original trait scale), \code{corrP},
#' \code{sigma_P} and the fixed effects \code{beta[trait, j]}. \code{sigma_P} are
#' the residual standard deviations, so \code{P}, \code{corrP} and \code{sigma_P}
#' describe the covariance of the residuals after the fixed effects, not the
#' marginal covariance of \code{Y}.
#'
#' Factor fits add quantities that describe the decomposition. \code{communality}
#' and \code{uniqueness} do not depend on how the factors are oriented, so they
#' are identified by the data. \code{loadings} and \code{factor_magnitude} do:
#' two solutions with the same \code{P} can differ by a rotation, a sign flip or
#' a permutation of the columns, and nothing fixes that choice except the prior.
#' A poor \code{Rhat} on those two is therefore not evidence that the fit
#' failed, and their posterior means describe whichever orientations the chains
#' visited. Judge convergence on \code{P}, \code{corrP}, \code{sigma_P},
#' \code{communality} and \code{uniqueness}, and pin the orientation down by
#' lowering \code{column_decay}, which orders the columns, or by asking for
#' fewer factors.
#' \itemize{
#'   \item \code{loadings}: the trait by factor loading matrix \eqn{\Lambda},
#'     on the correlation scale.
#'   \item \code{communality}: for trait \eqn{k},
#'     \eqn{c_k = \sum_m \lambda_{km}^2}, the part of the trait's variance
#'     explained by the common factors. On the correlation scale
#'     \eqn{c_k + u_k = 1}; multiplying by \code{sigma_P[k]^2} gives the common
#'     phenotypic variance of trait \eqn{k}.
#'   \item \code{uniqueness}: for trait \eqn{k}, \eqn{u_k = 1 - c_k}, the
#'     trait-specific variance left after the factors. Multiplying by
#'     \code{sigma_P[k]^2} gives the unique phenotypic variance.
#'   \item \code{factor_magnitude}: for factor \eqn{m},
#'     \eqn{\sum_k \lambda_{km}^2}, the variance the factor contributes across
#'     traits. It is a factor-level summary, not a per-trait quantity, and its
#'     values only compare across columns once the orientation is fixed.
#' }
#' These are available from \code{fit$summary()} and from
#' \code{\link{fit_draws}}: \code{loadings} is a trait by factor by draw array,
#' \code{factor_magnitude} is a factor by draw matrix, and the uniqueness and
#' communality vectors are trait by draw matrices.
#' @details The unstructured model estimates \code{P} as scales times an LKJ-distributed correlation matrix, with half-normal priors on the scales of the standardized responses. With \code{covariance_structure = "factor"}, \code{P} is instead \eqn{\mathrm{diag}(\sigma_P)(\Lambda\Lambda^T + \Psi)\mathrm{diag}(\sigma_P)}, with regularized-horseshoe loadings and deterministic column decay. The factor fit reports \code{loadings}, \code{factor_magnitude}, \code{uniqueness} and \code{communality} in addition to \code{P}, \code{corrP}, \code{sigma_P} and the fixed effects. Returned quantities are on the original response scale.
#'
#' Incomplete records are used rather than dropped. Without relatedness the
#' records are independent, so a row that misses some traits is evaluated on the
#' traits it does have: its marginal density is the multivariate normal on that
#' subvector, with the matching submatrix of \code{P}. Rows are grouped by which
#' traits they observe, each group gets its own Cholesky factor, and the missing
#' values are integrated out analytically, so they cost no parameters. A record
#' with nothing observed carries no likelihood and is an error instead. The count
#' of groups grows with the number of distinct patterns, at most one per row.
#'
#' The pointwise \code{log_lik} is the density the model block evaluates, on the
#' standardized response; the raw-scale density differs from it by a constant,
#' which does not affect comparisons between fits to the same data.
#'
#' The likelihood is computed by whitening the residuals with the Cholesky factor of \code{P} in a single triangular solve, at \code{O(N K^2 + K^3)} per evaluation rather than a per-individual loop, so the model scales to dozens of traits.
#'
#' Continuous predictors are centered and divided by their sample standard deviation inside Stan. Constant and 0/1 indicator columns are not centered or scaled. Standardized slopes have \code{normal(0, 1)} priors; the intercept contribution has a wider \code{normal(0, 2.5)} prior. Returned \code{beta} coefficients are backtransformed to original predictor and response units, including restoring the response mean to the intercept; \code{beta_scaled} draws are also available.
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 60
#' X <- cbind(Intercept = 1, sex = rbinom(n, 1, 0.5))
#' P <- random_covariance(3, min_var = 0.8, max_var = 1.2,
#'                        trait_names = c("t1", "t2", "t3"))
#' B <- rbind(c(10, 9, 11), c(1, -0.5, 0.3))
#' Y <- X %*% B + matrix(rnorm(n * 3), n) %*% chol(P)
#' colnames(Y) <- c("t1", "t2", "t3")
#'
#' fit <- lm_phen_cov(Y, X, iter = 1000, warmup = 500,
#'                    chains = 2, cores = 2, refresh = 0)
#' fit$summary(c("P", "corrP", "sigma_P"))
#' round(draws_mean(fit, "P"), 2)
#'
#' # Sparse factor P
#' factor_fit <- lm_phen_cov(Y, X, covariance_structure = "factor",
#'                           n_factors = 2, iter = 1000, warmup = 500,
#'                           chains = 2, cores = 2, refresh = 0)
#' factor_fit$summary(c("P", "loadings", "factor_magnitude"))
#' dim(fit_draws(factor_fit)$loadings)
#' }
#' @author Diogo Melo
#' @export
lm_phen_cov <- function(Y, X = NULL, lkj_prior = 2, iter = 3000, warmup = floor(iter/2),
                        chains = 4, cores = chains, control = list(adapt_delta = 0.99),
                        cmd_stan = TRUE, model_file = NULL,
                        covariance_structure = c("unstructured", "factor"),
                        n_factors = NULL, hs_global_scale = 1,
                        hs_slab_scale = 2, column_decay = 0.8, ...){
  lkj_prior_supplied <- !missing(lkj_prior)
  hs_global_scale_supplied <- !missing(hs_global_scale)
  hs_slab_scale_supplied <- !missing(hs_slab_scale)
  column_decay_supplied <- !missing(column_decay)
  covariance_structure <- match.arg(covariance_structure)
  factor_model <- covariance_structure == "factor"
  lkj_prior <- .check_positive(lkj_prior, "lkj_prior")
  hs_global_scale <- .check_positive(hs_global_scale, "hs_global_scale")
  hs_slab_scale <- .check_positive(hs_slab_scale, "hs_slab_scale")
  column_decay <- .check_in_range(column_decay, "column_decay", 0, 1,
                                  open_low = TRUE)
  cmd_stan <- .check_flag(cmd_stan, "cmd_stan")
  .check_iterations(iter, warmup)
  .check_sampler_size(chains, cores)
  YX <- .check_response_design(Y, X, allow_na = TRUE)
  Y <- YX$Y
  X <- YX$X
  missing <- YX$missing

  if (factor_model) {
    if (ncol(Y) < 2L)
      stop("`covariance_structure = \"factor\"` needs at least two traits.",
           call. = FALSE)
    if (!is.null(model_file))
      stop("`model_file` cannot be combined with `covariance_structure = \"factor\"`.",
           call. = FALSE)
    if (is.null(n_factors)) {
      n_factors <- .default_n_factors(ncol(Y))
      warning("`n_factors` was not supplied; using ", n_factors,
              " factors for ", ncol(Y), " traits. See `suggest_n_factors()` ",
              "for a number chosen from the data.", call. = FALSE)
    }
    n_factors <- .check_count(n_factors, "n_factors")
    if (n_factors >= ncol(Y))
      stop("`n_factors` must be smaller than the number of traits (",
           ncol(Y), ").", call. = FALSE)
  } else {
    if (!is.null(n_factors))
      stop("`n_factors` is only used with `covariance_structure = \"factor\"`.",
           call. = FALSE)
  }
  .warn_ignored(lkj_prior_supplied, !factor_model, "lkj_prior")
  .warn_ignored(hs_global_scale_supplied, factor_model, "hs_global_scale")
  .warn_ignored(hs_slab_scale_supplied, factor_model, "hs_slab_scale")
  .warn_ignored(column_decay_supplied, factor_model, "column_decay")

  cmd_stan <- .resolve_stan_backend(cmd_stan)

  stan_data = list(K = ncol(Y),
                   J = ncol(X),
                   N = nrow(Y),
                   X = X,
                   Y = Y,
                   miss = missing$miss,
                   n_pattern = missing$n_pattern,
                   pattern_obs = missing$pattern_obs,
                   # as.array: a one-element R vector is read as a scalar by
                   # rstan, and these are declared as arrays of the pattern count
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
  default_model <- if (factor_model) "phenotypicModelFactor.stan" else "phenotypicModel.stan"
  model_code <- .read_animal_model_code(model_file, default_model)

  if(cmd_stan){
    fit = cstan(model_code = model_code,
                data = stan_data,
                iter = iter, warmup = warmup, chains = chains, cores = cores,
                control = control, ...)
  } else {
    model_name <- if (factor_model) "phenotypicModelFactor" else "phenotypicModel"
    model <- .get_rstan_model(
      model_code,
      paste0(model_name, "_runtime"),
      if (is.null(model_file)) model_name else NULL
    )
    fit = .rstan_sample(model, stan_data, iter, warmup, chains, control, ...)
    fit <- .rename_fixed_effect_draws(fit, X, ncol(Y))
  }
  # the names are what fit_draws() labels its arrays with, so they have to be
  # recorded for both backends, not just the cmdstanr path
  attr(fit, "x_names") <- colnames(X)
  attr(fit, "y_names") <- colnames(Y)
  attr(fit, "covariance_structure") <- covariance_structure
  if (factor_model) attr(fit, "n_factors") <- n_factors
  attr(fit, "lkj_prior") <- lkj_prior
  attr(fit, "column_decay") <- column_decay
  attr(fit, "custom_model") <- !is.null(model_file)
  # what the predictive check needs: the response as passed in and the design
  attr(fit, "stanimal_data") <- list(family = "phenotypic", Y = YX$Y_raw,
                                     X = X, A = NULL)
  return(fit)
}
