# Wrapped so that the tests can pretend a package is not installed.
.is_installed <- function(package) requireNamespace(package, quietly = TRUE)

# Which Stan engine to use. Both cmdstanr and rstan are suggested packages, so
# a fit falls back to rstan when cmdstanr is not installed, and errors with the
# installation instructions when neither is available.
.resolve_stan_backend <- function(cmd_stan) {
  if (cmd_stan) {
    if (.is_installed("cmdstanr")) return(TRUE)
    if (!.is_installed("rstan"))
      stop("fitting a model needs cmdstanr or rstan: install one of them with ",
           "install.packages(\"cmdstanr\", repos = c(\"https://stan-dev.r-universe.dev\", ",
           "getOption(\"repos\"))) (then ",
           "cmdstanr::install_cmdstan()) or install.packages(\"rstan\").",
           call. = FALSE)
    warning("cmdstanr is not installed, falling back to rstan. Install it with ",
            "install.packages(\"cmdstanr\", repos = c(\"https://stan-dev.r-universe.dev\", ",
            "getOption(\"repos\"))), then run ",
            "cmdstanr::install_cmdstan().", call. = FALSE)
    return(FALSE)
  }
  if (!.is_installed("rstan"))
    stop("cmd_stan = FALSE needs rstan: install.packages(\"rstan\").",
         call. = FALSE)
  FALSE
}

#' Fit linear mixed animal model
#' 
#' Fits a univariate or multivariate linear mixed animal model using a known relatedness matrix
#' @param Y response variable, univariate or multivariate. Every trait must vary: constant traits cannot be standardized and are an error. Missing values (\code{NA}) are allowed and are not discarded: see Details.
#' @param X Optional finite numeric fixed-effects matrix with the same rows as Y and full column rank. The default \code{NULL} fits an intercept-only model: one column of ones, named \code{"(Intercept)"}. A matrix must be able to represent a constant (an intercept column, or the 0/1 indicators of every level of a factor), so that the mean of the response can be reported. Continuous columns are centered and scaled internally; returned beta draws use original X units.
#' @param A Finite symmetric positive-definite relationship matrix matching the rows of Y. A matrix that is only positive semidefinite (duplicated individuals, for example) is rejected; set \code{blend_A = TRUE} to regularize it. When \code{A} and \code{Y} both have row names, they are used to put \code{A} in the order of \code{Y}; names that disagree are an error. Without names the rows are matched by position.
#' @param lkj_prior Shape parameter of the LKJ prior on the correlation matrices of the multivariate model. \code{1} is uniform over correlation matrices, values above \code{1} shrink correlations toward zero, and values below \code{1} favor larger correlations. The default \code{2} is recommended for many traits: at high dimensionality the uniform prior concentrates mass near \code{+/-1} and overestimates the magnitude of genetic correlations, a well-known problem in quantitative genetics. A finite number above 0. It is used only by multivariate models with \code{residual_structure = "unstructured"} and not \code{shared_factor}; explicitly supplying it in other configurations warns that it is ignored.
#' @param beta_prior Parameters of the beta prior on h2. Used by the univariate model with \code{parametrization = "h2"}, by the unstructured multivariate model with \code{parametrization = "h2"}, by the ordinary factor model for the trait heritabilities, and by the shared-factor model for the factor heritabilities. The default \code{c(1, 1)} is flat on h2 and on each factor heritability. It is validated for every model; explicitly supplying it where it is not used warns that it is ignored.
#' @param parametrization Which parameters carry the prior. \code{"variance"} (default) estimates the genetic and residual standard deviations with half-normal priors and derives h2. \code{"h2"} estimates a phenotypic standard deviation and a heritability per trait, and puts the beta prior (\code{beta_prior}) on h2. It is available for the univariate models and for the unstructured multivariate model, and the two differ only in the prior they imply on h2: because the response is standardized trait by trait, \eqn{s_G^2 + s_E^2} is the phenotypic variance, so the half-normal priors of \code{"variance"} are priors on a circle and they induce \eqn{p(h^2) \propto 1/\sqrt{h^2(1-h^2)}}, a density whose mass sits at the two ends. The two parameterizations agree when the data are informative and separate when they are not, so prefer \code{"h2"} when there are few individuals relative to the number of traits, when h2 is the object of inference rather than \code{G}, or when a \code{"variance"} fit puts its heritabilities against 0 or 1. No sample-size threshold has been validated, so treat that as a reason to check rather than a rule. The ordinary factor model has no switch: its scales are locked to the phenotypic scale and it always puts the beta prior on the trait heritabilities. Supplying this argument for the factor or shared-factor model warns that it is ignored, and a univariate \code{"variance"} fit with low E-BFMI already warns and suggests \code{"h2"}.
#' @param iter Number of MCMC iterations, a whole number.
#' @param warmup Number of warmup iterations, a whole number smaller than \code{iter}.
#' @param chains Number of chains to run, a whole number.
#' @param cores Number of chains to run in parallel (cmdstanr only), a whole number.
#' @param model_file Optional path to a custom Stan file replacing the default model.
#' @param control Named list of sampler controls, passed to whichever engine runs: \code{adapt_delta}, \code{adapt_engaged}, \code{max_treedepth}, \code{step_size} (rstan's \code{stepsize} is accepted), \code{init_buffer}, \code{term_buffer}, \code{window} and \code{metric}. Any other name is an error, so a control is never silently ignored.
#' @param cmd_stan Logical, use cmdstanr. If FALSE, uses rstan. Both engines are suggested packages: when cmdstanr is not installed the fit falls back to rstan with a warning, and when neither is installed the fit errors.
#' @param blend_A If TRUE, replace A with 0.99 * A + 0.01 * I before fitting, making singular A positive definite.
#' @param residual_structure Residual covariance for multivariate fits: "unstructured" (default, full covariance between traits), "diagonal", or "shared_factor". With a diagonal residual, environmental covariances between traits are absorbed into the genetic covariance and inflate heritability, so "unstructured" is recommended when estimating G without a factor model. "shared_factor" requires \code{genetic_structure = "factor"} and fits the shared-factor model \eqn{G = \Lambda \mathrm{diag}(h^2_\text{factor}) \Lambda^T + \Psi_G}, \eqn{E = \Lambda \mathrm{diag}(1 - h^2_\text{factor}) \Lambda^T + \Psi_E}, so the same latent factors explain genetic and residual covariance and each factor has a heritability.
#' @param genetic_structure Genetic covariance for multivariate fits: "unstructured" (default, direct scales and LKJ correlation) or "factor", which uses the sparse factor model \eqn{G = \mathrm{diag}(s_G)(\Lambda\Lambda^T + \Psi)\mathrm{diag}(s_G)} on the correlation scale. The factor model gives loadings regularized-horseshoe priors for automatic shrinkage; loadings are not exactly zero. The ordinary factor model reports \code{loadings}, \code{factor_magnitude}, \code{uniqueness}, \code{communality} and \code{s_G}. With \code{residual_structure = "shared_factor"} the same loadings also structure E, and the fit reports \code{loadings}, \code{factor_magnitude}, \code{factor_heritability}, \code{genetic_uniqueness}, \code{residual_uniqueness}, \code{uniqueness} and \code{communality}. It is available only for multivariate fits and is not compatible with \code{model_file}.
#' @param n_factors Number of factors when \code{genetic_structure = "factor"}; for \code{residual_structure = "shared_factor"} these are shared by \code{G} and \code{E}. It is a whole number between 1 and the number of traits minus one. If \code{NULL}, the factor model falls back to \code{min(K - 1, max(2, floor(sqrt(K))), 20)} and warns with the value; that fallback depends only on the number of traits, so use \code{\link{suggest_n_factors}} to choose it from the data when the factors themselves are meant to be interpreted. Ignored otherwise.
#' @param hs_global_scale,hs_slab_scale Scales of the regularized-horseshoe prior on factor loadings: the global half-Cauchy scale and the slab scale. Smaller \code{hs_global_scale} gives stronger shrinkage and fewer effective loadings; \code{hs_slab_scale} bounds the largest loading in the slab. Both must be positive and are only used by \code{genetic_structure = "factor"}; explicitly supplying them elsewhere warns that they are ignored. See Piironen and Vehtari (2017) on choosing the global scale from the expected number of nonzero loadings.
#' @param column_decay Standard-deviation multiplier for successive factor columns, between 0 and 1. Column \code{m} is multiplied by \code{column_decay^(m - 1)}, so values below 1 favour decreasing factor magnitudes and make the factor labels interpretable. The default 0.8 orders the columns without forcing a strict hierarchy; with \code{column_decay = 1} the columns are exchangeable a priori and their labels are free to swap between chains. How strong the decay needs to be depends on how many columns have to be told apart: with two columns of comparable size 0.8 is often too mild and 0.5 is what separates them, while a strong decay works against a model with many columns, because \code{0.5^(m - 1)} cuts the prior of the fifth column to a sixteenth of the first. This is a deterministic prior on the column order, not a learned decay parameter; it is used by both factor structures, including the shared-factor model. Explicitly supplying it outside \code{genetic_structure = "factor"} warns that it is ignored.
#' @param ... additional parameters for stan sampling
#' @return With \code{cmd_stan = TRUE}, a \code{CmdStanMCMC} object from cmdstanr: use \code{fit$draws()} (a \code{posterior} draws array) and \code{fit$summary()}, and pass these, or the fit itself, to bayesplot (e.g. \code{bayesplot::mcmc_trace(fit$draws("h2"))}, \code{bayesplot::nuts_params(fit)}). Fixed effects are \code{beta[j]} (univariate) or \code{beta[trait,j]} for column \code{j} of \code{X}; use \code{\link{rename_beta}} to get draws named after the columns of \code{X}. Use \code{\link{fit_draws}} to reshape the draws into trait by trait by draw arrays, and \code{\link{draws_mean}} or \code{\link{draws_quantile}} for posterior summaries that are already matrices. Factor fits also expose \code{loadings}, \code{factor_magnitude}, \code{uniqueness} and \code{communality}; shared-factor fits add \code{factor_heritability}, \code{genetic_uniqueness} and \code{residual_uniqueness}. With \code{cmd_stan = FALSE}, a \code{stanfit} object whose \code{beta} draws are named after the columns of \code{X}.
#' Multivariate fits report \code{G}, \code{E}, \code{P = G + E}, \code{corrG},
#' \code{sigma_G}, \code{sigma_E} and \code{h2 = diag(G) / diag(P)}, plus the
#' fixed effects \code{beta[trait, j]}.
#'
#' Every fit reports \code{log_lik}, one term per row of the response, and the
#' terms sum to the joint log density. The basis is the catch: the model block
#' works on the eigen-rotation of \code{A}, where the records are independent, so
#' term \code{n} belongs to rotated record \code{n} and not to individual
#' \code{n}. With \code{A} equal to the identity the two coincide. Otherwise the
#' sum is still the exact joint log likelihood, but leave-one-out over rotated
#' records is not leave-one-individual-out, so treat \code{loo} results as a
#' model-comparison criterion rather than an individual-level predictive check.
#'
#' Factor fits add quantities that describe the factor decomposition. They come
#' in two groups. \code{communality}, \code{uniqueness} and
#' \code{factor_heritability} do not depend on how the factors are oriented, so
#' the data identify them and their \code{Rhat} and effective sample size are
#' read like any other parameter. \code{loadings} and \code{factor_magnitude}
#' do depend on it: two solutions with the same \code{G} can differ by a
#' rotation, a sign flip or a permutation of the columns, and nothing fixes that
#' choice except the prior. A poor \code{Rhat} on those two is therefore not
#' evidence that the fit failed, and their posterior means describe whichever
#' orientations the chains visited. Judge convergence on \code{G}, \code{E},
#' \code{P}, \code{corrG}, \code{communality} and \code{uniqueness}, and pin
#' the orientation down by lowering \code{column_decay}, which orders the
#' columns, or by asking for fewer factors.
#' \itemize{
#'   \item \code{loadings}: the trait by factor loading matrix \eqn{\Lambda}.
#'     For the ordinary factor model these are on the correlation scale; for
#'     shared-factor fits they are on the standardized phenotype scale.
#'   \item \code{communality}: for trait \eqn{k},
#'     \eqn{c_k = \sum_m \lambda_{km}^2}. This is the part of the trait's
#'     variance explained by the common factors. In the ordinary factor model
#'     it is on the correlation scale, so \eqn{c_k + u_k = 1}; multiplying by
#'     \eqn{s_{Gk}^2} gives the common genetic variance of trait \eqn{k}. In
#'     shared-factor fits it is the common phenotypic variance, and
#'     \eqn{c_k + u_k} equals the trait's standardized phenotypic variance.
#'   \item \code{uniqueness}: for trait \eqn{k}, the trait-specific variance
#'     left after the factors. In the ordinary factor model
#'     \eqn{u_k = 1 - c_k} is the unique genetic variance on the correlation
#'     scale; multiply by \eqn{s_{Gk}^2} for the unique genetic variance. In
#'     shared-factor fits \eqn{u_k = \psi_{Gk} + \psi_{Ek}}, with
#'     \code{genetic_uniqueness} \eqn{\psi_{Gk}} and
#'     \code{residual_uniqueness} \eqn{\psi_{Ek}} giving the two parts.
#'   \item \code{factor_magnitude}: for factor \eqn{m},
#'     \eqn{\sum_k \lambda_{km}^2}, the variance the factor contributes across
#'     traits. It is a factor-level summary, not a per-trait quantity. Ranking
#'     the columns only means something once the orientation is fixed, by
#'     \code{column_decay} or by an \code{n_factors} small enough to leave the
#'     factors identified.
#'   \item \code{factor_heritability}: shared-factor fits only. For factor
#'     \eqn{m}, \eqn{h^2_m} is the fraction of that factor's variance that is
#'     genetic: the genetic part is multiplied by \eqn{h^2_m} and the residual
#'     part by \eqn{1-h^2_m}.
#'   \item \code{s_G}: the genetic standard deviation of each trait on the
#'     standardized response scale, which is the diagonal scale that moves the
#'     correlation-scale decomposition to the original scale; multiplying it by
#'     the response standard deviation gives the genetic standard deviations
#'     \code{sigma_G} on the data scale. It is reported by every multivariate
#'     fit except the shared-factor model, which scales its factors by
#'     \code{sqrt_psi_G} and \code{sqrt_psi_E} instead. In the ordinary factor
#'     model it is derived from \code{h2} and \code{s_P}.
#' }
#' These are available from \code{fit$summary()} and from
#' \code{\link{fit_draws}}. The latter returns \code{loadings} as a trait by
#' factor by draw array, \code{factor_heritability} and
#' \code{factor_magnitude} as factor by draw matrices, and the uniqueness and
#' communality vectors as trait by draw matrices.
#'
#' @details By default the models estimate the additive genetic and residual standard deviations (\code{sigma_G}, \code{sigma_E}) directly, with half-normal priors on the standardized response scale, and report \code{h2} as a derived quantity. With \code{parametrization = "h2"}, and always in the ordinary factor model, the two scales are written as \eqn{s_G = \sqrt{h^2}s_P} and \eqn{s_E = \sqrt{1-h^2}s_P} instead, so \code{beta_prior} is a prior on the trait heritabilities and \eqn{s_P} is the phenotypic scale. The multivariate model estimates the additive genetic covariance matrix \code{G} and the residual covariance \code{E} directly, each as scales times an LKJ-distributed correlation matrix (\code{lkj_prior}), with half-normal priors on the scales of the standardized responses. With \code{genetic_structure = "factor"}, \code{G} is instead \eqn{\mathrm{diag}(s_G)(\Lambda\Lambda^T + \Psi)\mathrm{diag}(s_G)}, where the loadings \eqn{\Lambda} have regularized-horseshoe priors and \eqn{\Psi} is diagonal on the correlation scale; the default global shrinkage is \code{hs_global_scale = 1}. With \code{residual_structure = "shared_factor"}, the same loadings enter both covariances: \eqn{G = \Lambda\mathrm{diag}(h^2_\text{factor})\Lambda^T + \mathrm{diag}(\psi_G)} and \eqn{E = \Lambda\mathrm{diag}(1 - h^2_\text{factor})\Lambda^T + \mathrm{diag}(\psi_E)}, with the beta prior on the factor heritabilities controlled by \code{beta_prior}. The optional \code{column_decay} multiplies the prior standard deviation of factor column \eqn{m} by \code{column_decay^(m - 1)}, which favours factors in decreasing order without adding a learned decay parameter. Breeding values are marginalized using the eigendecomposition of \code{A}: the response and the design are rotated by its eigenvectors and the likelihood is computed on the rotated records, with the eigenvalues entering through the per-record covariance. No per-individual random effects are needed, which is why a singular \code{A} can be fitted once it is regularized with \code{blend_A = TRUE}. Returned quantities (original response scale) are \code{G}, \code{E}, \code{P = G + E}, \code{corrG}, \code{sigma_G}, \code{sigma_E} and \code{h2 = diag(G) / diag(P)}.
#'
#' With \code{parametrization = "h2"}, and always in the ordinary factor model, the scales are written as a phenotypic scale and a heritability rather than as two standard deviations, because the response is standardized trait by trait and the two scales are tied to its variance. That makes \eqn{s_G^2 + s_E^2} the phenotypic variance of the trait and fixes it near 1, so half-normal priors on \eqn{s_G} and on \eqn{s_E} are priors on a circle. The density they induce on \eqn{h^2} is proportional to \eqn{1/\sqrt{h^2(1-h^2)}}, an arcsine prior that puts about 14 percent of its mass below 0.05 and 14 percent above 0.95, which makes a moderate heritability the least likely answer a priori. With few individuals the likelihood is weak enough for that to show through, and a fit can then sit at either end of the interval. Writing the scales as \eqn{s_G = \sqrt{h^2}s_P} and \eqn{s_E = \sqrt{1-h^2}s_P} puts \code{beta_prior} on \eqn{h^2} directly, so the prior applied is the one that was asked for, and \eqn{s_P} is the well-determined quantity that is left.
#'
#' The per-individual likelihoods are whitened by the residual Cholesky factor and rotated into a single eigendecomposition of \code{E^{-1/2} G E^{-1/2}}, so each evaluation costs \code{O(N K^2 + K^3)} rather than \code{O(N K^3)}. This lets the multivariate model scale to dozens of traits without a Cholesky factorization per individual. When the eigenvalues of that matrix are nearly coincident (for example at initialization, or when \code{G} is roughly proportional to \code{E}) the model automatically falls back to a per-individual Cholesky evaluation, which has a well-behaved gradient.
#'
#' The cost that grows with the number of individuals is the dense relationship matrix: \code{A} is passed in full and each chain eigendecomposes it before sampling (\code{O(N^3)} time and \code{O(N^2)} memory per chain). That cost is paid once, so with many traits sampling dominates and the decomposition is a small share of the fit, while with few traits and short runs it can be most of it. The models are intended for many traits on hundreds to a few thousand individuals. For cohorts of tens of thousands, where that matrix is the bottleneck, an engine that uses the sparse inverse of \code{A} (such as ASReml) is the better tool.
#'
#' Incomplete records are used rather than dropped. A missing value is not
#' discarded, because with relatedness it carries information: the value is
#' informed by every relative's observed traits. The gaps are filled by
#' parameters on the standardized response scale, drawn with the rest of the
#' model and reported as \code{Y_mis}, and a flat prior on them makes the
#' posterior over the other parameters the exact marginal over the observed
#' data. Each gap is one parameter, so the cost grows with the number of gaps;
#' a record with nothing observed is allowed, since its relatives still inform
#' it. Dropping incomplete rows instead is the alternative, and it throws away
#' the rows that are complete for the other traits.
#'
#' Two configurations are worth avoiding when many values are missing. A
#' \code{residual_structure = "diagonal"} model with a small sample is weakly
#' identified once gaps loosen the residual structure, and can mix badly:
#' prefer \code{"unstructured"}, which is the default and what the missing-data
#' machinery is tested against. And \code{loo} on a fit with gaps compares
#' models over the rotated records, as it does for complete data, so the terms
#' are not leave-one-individual-out.
#'
#' The LKJ shape \code{lkj_prior} (default 2) deliberately shrinks correlations toward zero. With many traits, an unstructured correlation matrix fit under a flat prior tends to overestimate the magnitude of correlations, especially genetic correlations between traits with limited information, so a shape above 1 acts as a mild regularizer; increase it for stronger shrinkage. Note that fitting many traits remains expensive: the likelihood is cheap, but an unstructured \code{G} for \code{K = 40} traits has on the order of 1700 parameters and NUTS needs many leapfrog steps, so a single fit can take an hour or more. Fewer iterations, a larger \code{lkj_prior}, or \code{residual_structure = "diagonal"} all make fits faster.
#'
#' Continuous predictors are centered and divided by their sample standard deviation inside Stan. Constant and 0/1 indicator columns are not centered or scaled. Standardized slopes have \code{normal(0, 1)} priors; the intercept contribution has a wider \code{normal(0, 2.5)} prior. Returned \code{beta} coefficients are backtransformed to original predictor and response units, including restoring the response mean to the intercept; \code{beta_scaled} draws are also available.
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 40
#' ped <- data.frame(
#'   id = 1:n,
#'   sire = c(rep(NA, 20), sample(1:10, 20, replace = TRUE)),
#'   dam = c(rep(NA, 20), sample(11:20, 20, replace = TRUE))
#' )
#' A <- makeA(ped)
#' X <- cbind(Intercept = 1, sex = rbinom(n, 1, 0.5))
#' G <- random_covariance(2, min_var = 0.8, max_var = 1.2,
#'                        trait_names = c("t1", "t2"))
#' E <- random_covariance(2, min_var = 0.8, max_var = 1.2,
#'                        trait_names = c("t1", "t2"))
#' Y <- X %*% rbind(c(10, 9), c(1, -0.5)) +
#'   t(chol(A)) %*% matrix(rnorm(n * 2), n) %*% chol(G) +
#'   matrix(rnorm(n * 2), n) %*% chol(E)
#' colnames(Y) <- c("t1", "t2")
#'
#' fit <- lmm_animal(Y, X, A, iter = 1000, warmup = 500,
#'                   chains = 2, cores = 2, refresh = 0)
#' fit$summary(c("G", "E", "h2"))
#' round(draws_mean(fit, "G"), 2)
#'
#' # Univariate model with the h2 parametrization
#' fit_h2 <- lmm_animal(Y[, 1], X, A, parametrization = "h2",
#'                      iter = 1000, warmup = 500,
#'                      chains = 2, cores = 2)
#'
#' # Sparse factor G
#' factor_fit <- lmm_animal(Y, X, A, genetic_structure = "factor",
#'                          n_factors = 1, iter = 1000, warmup = 500,
#'                          chains = 2, cores = 2)
#' factor_fit$summary(c("loadings", "factor_magnitude"))
#' dim(fit_draws(factor_fit)$loadings)
#'
#' # Shared factors for G and E
#' shared_fit <- lmm_animal(Y, X, A, genetic_structure = "factor",
#'                          residual_structure = "shared_factor",
#'                          n_factors = 1, iter = 1000, warmup = 500,
#'                          chains = 2, cores = 2)
#' shared_fit$summary(c("factor_heritability", "genetic_uniqueness"))
#' }
#' @export
#' @author Diogo Melo
lmm_animal = function(Y, X = NULL, A, lkj_prior = 2, beta_prior = c(1, 1),
                      parametrization = c("variance", "h2"),
                      iter = 3000, warmup = floor(iter/2), 
                      chains = 4, cores = chains, model_file = NULL, control = list(adapt_delta = 0.99), 
                      cmd_stan = TRUE, blend_A = FALSE,
                      residual_structure = c("unstructured", "diagonal",
                                             "shared_factor"),
                      genetic_structure = c("unstructured", "factor"),
                      n_factors = NULL, hs_global_scale = 1,
                      hs_slab_scale = 2, column_decay = 0.8, ...){
  beta_prior_supplied <- !missing(beta_prior)
  lkj_prior_supplied <- !missing(lkj_prior)
  hs_global_scale_supplied <- !missing(hs_global_scale)
  hs_slab_scale_supplied <- !missing(hs_slab_scale)
  column_decay_supplied <- !missing(column_decay)
  parametrization_supplied <- !missing(parametrization)
  residual_structure <- match.arg(residual_structure)
  genetic_structure <- match.arg(genetic_structure)
  parametrization <- match.arg(parametrization)
  lkj_prior <- .check_positive(lkj_prior, "lkj_prior")
  beta_prior <- .check_beta_prior(beta_prior)
  hs_global_scale <- .check_positive(hs_global_scale, "hs_global_scale")
  hs_slab_scale <- .check_positive(hs_slab_scale, "hs_slab_scale")
  column_decay <- .check_in_range(column_decay, "column_decay", 0, 1,
                                  open_low = TRUE)
  blend_A <- .check_flag(blend_A, "blend_A")
  cmd_stan <- .check_flag(cmd_stan, "cmd_stan")
  .check_iterations(iter, warmup)
  .check_sampler_size(chains, cores)
  YX <- .check_response_design(Y, X, allow_na = TRUE, allow_empty_rows = TRUE)
  Y <- YX$Y
  X <- YX$X
  missing <- YX$missing
  A <- as.matrix(A)

  factor_model <- genetic_structure == "factor"
  if (residual_structure == "shared_factor" && !factor_model)
    stop("`residual_structure = \"shared_factor\"` requires `genetic_structure = \"factor\"`.",
         call. = FALSE)
  shared_factor_model <- factor_model && residual_structure == "shared_factor"
  # The ordinary factor model has no choice about its parameters, so this only
  # distinguishes the two parameterisations of the unstructured multivariate
  # model. Trait h2 is what the beta prior is about either way.
  h2_param <- ncol(Y) > 1L && !factor_model && parametrization == "h2"
  if (factor_model) {
    if (ncol(Y) < 2L)
      stop("`genetic_structure = \"factor\"` needs at least two traits.",
           call. = FALSE)
    if (!is.null(model_file))
      stop("`model_file` cannot be combined with `genetic_structure = \"factor\"`.",
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
      stop("`n_factors` is only used with `genetic_structure = \"factor\"`.",
           call. = FALSE)
  }

  beta_used <- (ncol(Y) == 1L && parametrization == "h2") || factor_model ||
    h2_param
  lkj_used <- ncol(Y) > 1L && !shared_factor_model &&
    residual_structure == "unstructured"
  factor_used <- factor_model
  .warn_ignored(beta_prior_supplied, beta_used, "beta_prior")
  .warn_ignored(lkj_prior_supplied, lkj_used, "lkj_prior")
  # the factor model has no choice of parameterisation, so the argument applies
  # to every other fit, univariate or not
  .warn_ignored(parametrization_supplied, !factor_model, "parametrization")
  .warn_ignored(hs_global_scale_supplied, factor_used, "hs_global_scale")
  .warn_ignored(hs_slab_scale_supplied, factor_used, "hs_slab_scale")
  .warn_ignored(column_decay_supplied, factor_used, "column_decay")

  n <- nrow(Y)
  A <- .check_relationship_matrix(A, Y)
  A <- .check_symmetric(A)
  eigenvalues <- .check_positive_semidefinite(A)

  if (blend_A) {
    A <- 0.99 * A + 0.01 * diag(n)
    eigenvalues <- .relationship_eigenvalues(A)
  }
  if (min(eigenvalues) <= .eigen_tolerance(eigenvalues))
    stop("`A` must be positive definite for Stan; set `blend_A = TRUE` to ",
         "regularize it.", call. = FALSE)

  cmd_stan <- .resolve_stan_backend(cmd_stan)
      
  if(ncol(Y) == 1L){
    stan_data = list(J = ncol(X),
                     N = nrow(Y),
                     X = X,
                     Y = as.numeric(Y),
                     A = A,
                     miss = as.integer(missing$miss),
                     n_missing = missing$n_missing,
                     mis_row = as.array(missing$mis_row),
                     x_mean_weight = YX$x_mean_weight)
    uni_file <- if(parametrization == "h2") "animalModelUniH2.stan" else "animalModelUni.stan"
    if(parametrization == "h2")
      stan_data$beta_prior <- beta_prior
    model_code <- .read_animal_model_code(model_file, uni_file)
    if(cmd_stan){
      fit = cstan(model_code = model_code,
                  data = stan_data,
                  iter = iter, warmup = warmup, chains = chains, cores = cores,
                  control = control, ...)
    }
    else {
      model <- .get_rstan_model(
        model_code,
        paste0(sub("\\.stan$", "", uni_file), "_runtime"),
        if (is.null(model_file)) sub("\\.stan$", "", uni_file) else NULL
      )
      fit = .rstan_sample(model, stan_data, iter, warmup, chains, control, ...)
      fit <- .rename_fixed_effect_draws(fit, X, ncol(Y))
    }
    if(parametrization == "variance") {
      .warn_low_bfmi(fit)
      .warn_h2_boundary(fit)
    }
  } else{
    stan_data = list(K = ncol(Y),
                     J = ncol(X),
                     N = nrow(Y),
                     X = X,
                     Y = Y,
                     A = A,
                     miss = missing$miss,
                     n_missing = missing$n_missing,
                     # as.array: a one-element R vector is read as a scalar by
                     # rstan, and these are declared as arrays of the gap count
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
      stan_data$unstructured_E <- as.integer(residual_structure == "unstructured")
      if (h2_param) stan_data$beta_prior <- beta_prior
      if (factor_model) {
        stan_data$M <- n_factors
        stan_data$hs_global_scale <- hs_global_scale
        stan_data$hs_slab_scale <- hs_slab_scale
        stan_data$column_decay <- column_decay
        # same name and meaning as in animalModelUniH2.stan: the beta prior on
        # the trait heritabilities. The shared-factor model above calls the same
        # pair h2_prior_a/b because there it is a factor heritability.
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
    if(cmd_stan){
      fit = cstan(model_code = model_code,
                  data = stan_data,
                  iter = iter, warmup = warmup, chains = chains , cores = cores,
                  control=control , ...)
    } else {
      model_name <- if (shared_factor_model) {
        "animalModelGSharedFactor"
      } else if (factor_model) {
        "animalModelGFactor"
      } else if (h2_param) {
        "animalModelGH2"
      } else {
        "animalModelG"
      }
      model <- .get_rstan_model(
        model_code,
        paste0(model_name, "_runtime"),
        if (is.null(model_file)) model_name else NULL
      )
      fit = .rstan_sample(model, stan_data, iter, warmup, chains, control, ...)
      fit <- .rename_fixed_effect_draws(fit, X, ncol(Y))
    }
    # Only the variance parametrization of the unstructured model implies the
    # arcsine prior this check is about. The ordinary factor model already puts a
    # beta prior on h2, and the shared-factor model puts one on the factor
    # heritabilities, so there is nothing to suggest for either.
    if (!factor_model && parametrization == "variance")
      .warn_h2_boundary(fit)
  }
  # the names are what fit_draws() labels its arrays with, so they have to be
  # recorded for both backends, not just the cmdstanr path
  attr(fit, "x_names") <- colnames(X)
  attr(fit, "y_names") <- colnames(Y)
  attr(fit, "genetic_structure") <- genetic_structure
  attr(fit, "residual_structure") <- residual_structure
  # the parametrization only means something where the model has a choice, so
  # record what the model actually used rather than what was passed. The
  # ordinary factor model has none: its two scales are tied to the phenotypic
  # scale, so it always puts the prior on h2, and "h2" is what it ran.
  attr(fit, "parametrization") <- if (ncol(Y) == 1L) {
    parametrization
  } else if (factor_model && !shared_factor_model) {
    "h2"
  } else if (h2_param) {
    "h2"
  } else {
    "variance"
  }
  attr(fit, "lkj_prior") <- lkj_prior
  attr(fit, "beta_prior") <- beta_prior
  # a custom Stan file can use any prior it likes, so nothing below can be
  # recovered from the recorded arguments alone
  attr(fit, "custom_model") <- !is.null(model_file)
  # record every prior the fit ran with, used or not, so that the two fitters
  # expose the same metadata and a stored fit can be reconstructed without the
  # call that produced it
  attr(fit, "hs_global_scale") <- hs_global_scale
  attr(fit, "hs_slab_scale") <- hs_slab_scale
  attr(fit, "column_decay") <- column_decay
  if (factor_model) attr(fit, "n_factors") <- n_factors
  # what the predictive check needs, and what a breeding-value pass will need
  # later: the response as passed in, the design, and the relationship matrix
  attr(fit, "stanimal_data") <- list(family = "animal", Y = YX$Y_raw, X = X,
                                     A = A)
  return(fit)
}

# Shared by every model function: Y and X, and the weights that restore the
# response mean through the design. The Gibbs sampler does not shift the mean
# through X, so it passes require_constant = FALSE.
#
# X is optional. Without one the model is intercept only: a single column of
# ones, named as model.matrix() names the intercept, so the beta draws and the
# reported means come out of the same code path as an explicit intercept column.
#
# allow_na is set by the model functions that can handle incomplete traits, and
# those use the `missing` entry of the result to tell Stan where the gaps are.
# allow_empty_rows goes with it: the animal models impute from relatives, so a
# record with nothing observed is still informative.
.check_response_design <- function(Y, X, require_constant = TRUE,
                                   allow_na = FALSE, allow_empty_rows = FALSE) {
  Y <- .check_numeric_matrix(Y, "Y", allow_na = allow_na)
  Y_raw <- Y # with the gaps still in place, for the predictive checks
  if (is.null(X))
    X <- matrix(1, nrow(Y), 1L, dimnames = list(NULL, "(Intercept)"))
  X <- .check_numeric_matrix(X, "X")
  p <- ncol(X)
  if (nrow(Y) < 2L || p < 1L || nrow(X) != nrow(Y))
    stop("`Y` and `X` must have the same number of rows, and `X` must have ",
         "columns.", call. = FALSE)
  if (qr(X)$rank < p)
    stop("`X` must have full column rank.", call. = FALSE)
  missing <- .missing_response(Y, allow_empty_rows = allow_empty_rows)
  Y <- missing$Y
  .check_traits_vary(Y, missing$miss)
  weights <- NULL
  if (require_constant) {
    weights <- .intercept_weights(X)
    if (is.null(weights))
      stop("`X` cannot represent a constant, so the mean of the response ",
           "cannot be reported: include an intercept column, or the 0/1 ",
           "indicator columns of every level of a factor.", call. = FALSE)
  }
  list(Y = Y, X = X, x_mean_weight = weights, missing = missing,
       Y_raw = Y_raw)
}

# A trait without variance cannot be standardized, and gives no variance to
# split into genetic and residual parts. `miss` is the indicator matrix from
# .missing_response(), so each trait is judged on its observed values only.
.check_traits_vary <- function(Y, miss = NULL) {
  observed_sd <- if (is.null(miss)) {
    apply(Y, 2L, stats::sd)
  } else {
    vapply(seq_len(ncol(Y)),
           function(k) {
             present <- miss[, k] == 0L
             if (sum(present) < 2L) 0 else stats::sd(Y[present, k])
           },
           numeric(1))
  }
  constant <- observed_sd <= 0
  if (!any(constant))
    return(invisible(NULL))
  stop("every trait in `Y` must vary: ",
       paste(.trait_labels(Y)[constant], collapse = ", "),
       ". A constant trait has no additive or residual variance to estimate.",
       call. = FALSE)
}

# A symmetric relationship matrix, symmetrized against rounding.
.check_symmetric <- function(A, arg = "A") {
  scale <- max(1, max(abs(A)))
  if (max(abs(A - t(A))) > 1e-8 * scale)
    stop("`", arg, "` must be symmetric.", call. = FALSE)
  (A + t(A)) / 2
}

.relationship_eigenvalues <- function(A) {
  eigen(A, symmetric = TRUE, only.values = TRUE)$values
}

.eigen_tolerance <- function(values) 1e-8 * max(1, max(abs(values)))

# Returns the eigenvalues, which the callers reuse for the positive definite
# check after blending and, in the Gibbs sampler, to whiten the breeding values.
.check_positive_semidefinite <- function(A, arg = "A") {
  values <- .relationship_eigenvalues(A)
  if (min(values) < -.eigen_tolerance(values))
    stop("`", arg, "` must be positive semidefinite.", call. = FALSE)
  values
}

# Weights w with X w equal to a column of ones, when the constant is in the
# column space of X, and NULL when it is not. The models centre the responses,
# so the response mean is reported back through these weights: an intercept
# column gets all of it, the indicators of a factor share it.
.intercept_weights <- function(X) {
  weights <- qr.solve(X, rep(1, nrow(X)))
  if (anyNA(weights)) return(NULL)
  fitted <- drop(X %*% weights)
  if (max(abs(fitted - 1)) > 1e-8 * max(1, max(abs(fitted)))) return(NULL)
  weights
}

# Checks that A matches the rows of Y and returns it with its individuals in the
# order of Y. When both A and Y are named the names are used to align A; names
# that disagree are an error. Without names the rows are matched by position, so
# only the dimensions are checked.
.check_relationship_matrix <- function(A, Y) {
  Y <- as.matrix(Y)
  n <- nrow(Y)
  if (!is.numeric(A) || any(!is.finite(A)) || is.null(dim(A)))
    stop("A must be a finite numeric square matrix matching the rows of Y.",
         call. = FALSE)
  if (!identical(dim(A), c(n, n)))
    stop("A must be a square matrix matching the rows of Y: A is ", nrow(A),
         " x ", ncol(A), " but Y has ", n, " rows.", call. = FALSE)

  y_ids <- rownames(Y)
  a_ids <- rownames(A)
  if (is.null(y_ids) || is.null(a_ids))
    return(A)
  if (is.null(colnames(A)) || !identical(a_ids, colnames(A)))
    stop("A must have identical row and column names.", call. = FALSE)
  if (anyDuplicated(a_ids))
    stop("A has duplicated individual names.", call. = FALSE)

  missing_in_A <- setdiff(y_ids, a_ids)
  missing_in_Y <- setdiff(a_ids, y_ids)
  if (length(missing_in_A) || length(missing_in_Y))
    stop("the individuals in A do not match the rows of Y: ",
         length(missing_in_A), " rows of Y are not in A and ",
         length(missing_in_Y), " individuals in A are not in Y. Rename or ",
         "reorder them so the names agree, or drop the names from both to ",
         "match by position.", call. = FALSE)

  if (identical(y_ids, a_ids))
    return(A)
  A[y_ids, y_ids, drop = FALSE]
}

.warn_low_bfmi <- function(fit, threshold = 0.3) {
  bfmi <- tryCatch(
    if (inherits(fit, "stanfit")) rstan::get_bfmi(fit) else
      fit$diagnostic_summary(quiet = TRUE)$ebfmi,
    error = function(e) NULL)
  if (length(bfmi) && any(bfmi < threshold, na.rm = TRUE))
    warning(sum(bfmi < threshold, na.rm = TRUE), " of ", length(bfmi),
            " chains had E-BFMI below ", threshold, ". The variance ",
            "parametrization can sample poorly when the genetic variance is ",
            "weakly identified; consider parametrization = \"h2\".",
            call. = FALSE)
  invisible(fit)
}

# The variance parametrization puts half-normal priors on two standard
# deviations that are constrained to the phenotypic scale, which is a prior on a
# circle and induces p(h2) proportional to 1 / sqrt(h2 (1 - h2)) on each trait.
# That density is about three times higher at h2 = 0.05 than at 0.5, so a
# posterior mean against either boundary is compatible with the prior having
# decided the answer. The check reports the count rather than a verdict, because
# a trait with no additive variance legitimately sits at zero.
.warn_h2_boundary <- function(fit, low = 0.05, high = 0.95) {
  h2 <- tryCatch({
    d <- if (inherits(fit, "stanfit"))
      posterior::as_draws_matrix(fit) else
        fit$draws("h2", format = "draws_matrix")
    d[, grep("^h2\\[", posterior::variables(d)), drop = FALSE]
  }, error = function(e) NULL)
  if (is.null(h2) || !ncol(h2))
    return(invisible(fit))
  mean_h2 <- colMeans(as.matrix(h2))
  at_edge <- mean_h2 < low | mean_h2 > high
  if (any(at_edge))
    warning("h2 has a posterior mean below ", low, " or above ", high, " for ",
            sum(at_edge), " of ", length(mean_h2), " traits. The variance ",
            "parametrization implies a Beta(1/2, 1/2) prior on h2, which is ",
            "densest near 0 and near 1, so an answer at a boundary may be the ",
            "prior rather than the data: refit with parametrization = \"h2\" ",
            "before reading the heritabilities.", call. = FALSE)
  invisible(fit)
}

.beta_names <- function(x_names, n_traits) {
  if (n_traits == 1L)
    return(list(old = paste0("beta[", seq_along(x_names), "]"), new = paste0("beta[", x_names, "]")))
  grid <- expand.grid(trait = seq_len(n_traits), j = seq_along(x_names))
  list(old = paste0("beta[", grid$trait, ",", grid$j, "]"),
       new = paste0("beta[", grid$trait, ",", x_names[grid$j], "]"))
}

.rename_fixed_effect_draws <- function(fit, X, n_traits) {
  if (is.null(colnames(X)))
    return(fit)

  beta_indices <- grep("^beta\\[", names(fit))
  beta_names <- .beta_names(colnames(X), n_traits)$new
  if (length(beta_indices) == length(beta_names))
    names(fit)[beta_indices] <- beta_names
  fit
}

#' Posterior draws with the fixed effects named
#'
#' Extracts draws from a \code{lmm_animal} fit made with \code{cmd_stan = TRUE},
#' renaming the \code{beta} variables after the columns of \code{X}: named
#' \code{beta[name]} for one trait and \code{beta[trait,name]} for several,
#' instead of \code{beta[j]} and \code{beta[k,j]}. Every other variable is
#' returned unchanged, and if \code{X} had no column names the draws are
#' returned as \code{fit$draws()} gives them.
#' @param fit A \code{CmdStanMCMC} object returned by \code{lmm_animal}.
#' @param variables Optional character vector of variables to extract, using
#'   Stan's names (e.g. \code{"beta"}, \code{"h2"}). Default is all variables.
#' @param ... Further arguments passed to \code{fit$draws()}, e.g. \code{format}.
#' @return A \code{posterior} draws object, usable with bayesplot.
#' @export
rename_beta <- function(fit, variables = NULL, ...) {
  draws <- fit$draws(variables = variables, ...)
  x_names <- attr(fit, "x_names")
  if (is.null(x_names))
    return(draws)
  vars <- posterior::variables(draws)
  n_traits <- if (any(grepl("^beta\\[[0-9]+,[0-9]+\\]$", vars))) {
    max(as.integer(sub("^beta\\[([0-9]+),.*", "\\1", grep("^beta\\[[0-9]+,", vars, value = TRUE))))
  } else 1L
  map <- .beta_names(x_names, n_traits)
  if (n_traits == 1L && any(grepl("^beta\\[1,[0-9]+\\]$", vars)))
    map$old <- paste0("beta[1,", seq_along(x_names), "]")
  hit <- match(vars, map$old)
  vars[!is.na(hit)] <- map$new[hit[!is.na(hit)]]
  posterior::variables(draws) <- vars
  draws
}

.read_animal_model_code <- function(model_file, default = "animalModelG.stan") {
  if (!is.null(model_file))
    return(readChar(model_file, file.info(model_file)$size))

  model_filename <- default
  model_path <- system.file("stan", model_filename, package = "stanimal")
  if (!nzchar(model_path)) {
    candidates <- file.path(c("inst/stan", "package/inst/stan"), model_filename)
    model_path <- candidates[file.exists(candidates)][1L]
  }
  if (is.na(model_path) || !file.exists(model_path))
    stop("Could not locate the requested Stan model source.", call. = FALSE)
  paste(readLines(model_path, warn = FALSE), collapse = "\n")
}

.rstan_model_cache <- new.env(parent = emptyenv())

.get_rstan_model <- function(model_code, model_name, cache_key = NULL) {
  if (!is.null(cache_key) && exists(cache_key, envir = .rstan_model_cache,
                                    inherits = FALSE)) {
    cached <- get(cache_key, envir = .rstan_model_cache, inherits = FALSE)
    if (identical(cached$code, model_code))
      return(cached$model)
    model_name <- basename(tempfile(paste0(model_name, "_")))
  } else if (is.null(cache_key)) {
    model_name <- basename(tempfile(paste0(model_name, "_")))
  }

  model <- rstan::stan_model(model_code = model_code, model_name = model_name)
  if (!is.null(cache_key))
    assign(cache_key, list(code = model_code, model = model),
           envir = .rstan_model_cache)
  model
}

.rstan_sample <- function(model, stan_data, iter, warmup, chains, control, ...) {
  # Each chain draws its own start, as in cstan(): a shared start would make the
  # chains agree by construction. Overridden by an explicit `init` in `...`.
  control <- .rstan_control(.check_control(control))
  args <- list(object = model, data = stan_data, iter = iter, warmup = warmup,
               chains = chains, control = control)
  do.call(rstan::sampling, utils::modifyList(args, list(...)))
}
