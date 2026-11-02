#' Eigenvalues of the residual trait correlation matrix
#'
#' The eigenvalues behind a scree plot, for choosing \code{n_factors} in the
#' sparse factor models. The traits are residualized on \code{X} first, so the
#' eigenvalues describe the correlation the factor model actually fits.
#'
#' @param Y Numeric matrix of traits, individuals in rows.
#' @param X Optional finite numeric fixed-effects matrix with the same rows as
#'   \code{Y}, as in \code{\link{lmm_animal}}. The default \code{NULL} uses an
#'   intercept only, which is the same as centering the traits.
#' @return A data frame of class \code{scree} with one row per component and
#'   columns \code{component}, \code{eigenvalue}, \code{share} and
#'   \code{cumulative}. \code{share} is the eigenvalue divided by the number of
#'   traits and \code{cumulative} is its running total.
#' @details A factor contributes at most one eigenvalue to this decomposition,
#'   so the number of eigenvalues that stand out above the noise floor is an
#'   upper bound on a useful \code{n_factors}. Use
#'   \code{\link{suggest_n_factors}} to turn that into a number, and
#'   \code{\link{plot.scree}} to look at it.
#' @seealso \code{\link{suggest_n_factors}}, \code{\link{pc1_variance_share}}
#' @examples
#' set.seed(1)
#' P <- random_covariance(6, method = "factor", integration = 0.8)
#' Y <- matrix(stats::rnorm(300 * 6), 300) %*% chol(P)
#' colnames(Y) <- paste0("t", 1:6)
#' scree <- scree_eigenvalues(Y)
#' scree
#' plot(scree)
#' @export
scree_eigenvalues <- function(Y, X = NULL) {
  YX <- .check_response_design(Y, X, require_constant = FALSE)
  values <- .residual_eigenvalues(YX$Y, YX$X, traits = colnames(YX$Y))
  .scree_frame(values)
}

#' Number of factors suggested by parallel analysis
#'
#' Compares the eigenvalues of the residual trait correlation matrix with those
#' of data that has no correlation at all, and returns how many components stand
#' above that reference. It is the rule of thumb for setting \code{n_factors}
#' when you intend to interpret the factors.
#' @inheritParams scree_eigenvalues
#' @param nsim Number of simulated data sets, a whole number. More is slower
#'   and steadier; the default is enough to place the 95\% point.
#' @param level Quantile of the simulated eigenvalues to compare against, one
#'   number between 0 and 1. The default 0.95 keeps a component only when it
#'   clears 95\% of what uncorrelated data of the same size produces.
#' @param max_factors Largest number to report, a whole number. It is capped at
#'   the number of traits minus one, because the factor models need at least one
#'   residual degree of freedom per trait.
#' @return A whole number between 1 and \code{ncol(Y) - 1}. The raw count is in
#'   \code{attr(x, "retained")}, which is 0 when nothing stands above the noise,
#'   the eigenvalue cutoffs are in \code{attr(x, "reference")} for
#'   \code{\link{plot.scree}}, and \code{attr(x, "share")} is the fraction of
#'   the variance taken by the retained components.
#' @details Parallel analysis (Horn 1965) asks what the eigenvalues would look
#'   like with no shared structure: it simulates independent normal traits with
#'   the same number of individuals, removes the same fixed effects, and records
#'   the \code{level} quantile of each eigenvalue. A component is retained when
#'   its observed eigenvalue is larger. Unlike "keep the eigenvalues above one"
#'   or eyeballing the elbow, the cutoff then accounts for the sample size, the
#'   number of traits and the fixed effects.
#'
#'   The result is an \strong{upper bound}, not a recommendation to use that many
#'   factors. With a few hundred individuals and dozens of traits the noise floor
#'   is high enough that weak components clear it: on the mouse morphometric data
#'   it keeps seven dimensions, of which the first two take 28\% and 11\% of the
#'   variance while the seventh clears its cutoff by 0.003. Factors you intend to
#'   interpret are usually fewer: look at \code{\link{scree_eigenvalues}} and stop
#'   where the eigenvalues flatten out. The count also cannot tell a single strong
#'   direction split across two columns from two directions, which is what
#'   \code{column_decay} is for.
#'
#'   The result is at least 1 because the factor models need at least one factor;
#'   when nothing clears the reference, that is a sign that a factor model is
#'   fitting noise rather than structure, and \code{attr(x, "retained")} says so.
#'
#'   Fix the seed before calling when the answer has to be reproducible.
#' @seealso \code{\link{scree_eigenvalues}}, \code{\link{plot.scree}}
#' @examples
#' set.seed(1)
#' # two modules of three traits each
#' loadings <- matrix(0, 6, 2)
#' loadings[1:3, 1] <- 0.8
#' loadings[4:6, 2] <- 0.8
#' P <- loadings %*% t(loadings) + diag(6) * 0.36
#' Y <- matrix(stats::rnorm(400 * 6), 400) %*% chol(P)
#' colnames(Y) <- paste0("t", 1:6)
#' suggest_n_factors(Y)                     # 2
#' scree <- scree_eigenvalues(Y)
#' plot(scree, reference = suggest_n_factors(Y))
#' @export
suggest_n_factors <- function(Y, X = NULL, nsim = 200, level = 0.95,
                              max_factors = 20) {
  nsim <- .check_count(nsim, "nsim")
  max_factors <- .check_count(max_factors, "max_factors")
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1)
    stop("`level` must be one number between 0 and 1.", call. = FALSE)
  YX <- .check_response_design(Y, X, require_constant = FALSE)
  Y <- YX$Y
  X <- YX$X
  if (ncol(Y) < 2L)
    stop("`Y` must have at least two traits: a single trait leaves no factor "
         , "for the model to fit.", call. = FALSE)
  traits <- colnames(Y)
  observed <- .residual_eigenvalues(Y, X, traits = traits)
  p <- length(observed)
  n <- nrow(Y)
  simulated <- vapply(seq_len(nsim), function(i) {
    Z <- matrix(stats::rnorm(n * p), n, p)
    .residual_eigenvalues(Z, X, traits = traits)
  }, numeric(p))
  reference <- apply(simulated, 1L, stats::quantile, probs = level,
                     names = FALSE)
  retained <- sum(observed > reference)
  out <- as.integer(max(1L, min(retained, max_factors, p - 1L)))
  attr(out, "retained") <- as.integer(retained)
  attr(out, "eigenvalues") <- observed
  attr(out, "reference") <- reference
  attr(out, "level") <- level
  attr(out, "share") <- sum(observed[seq_len(out)]) / sum(observed)
  out
}

#' Scree plot
#'
#' @param x A \code{scree} object from \code{\link{scree_eigenvalues}}.
#' @param reference Optional eigenvalue cutoffs, either the object returned by
#'   \code{\link{suggest_n_factors}} or its \code{reference} attribute.
#'   Components above the dashed line are the ones parallel analysis keeps.
#' @param title Plot title.
#' @param ... Unused, for consistency with \code{plot}.
#' @return A ggplot object.
#' @seealso \code{\link{scree_eigenvalues}}, \code{\link{suggest_n_factors}}
#' @export
plot.scree <- function(x, reference = NULL, title = "Scree plot", ...) {
  # the object suggest_n_factors() returns carries the cutoffs as an attribute,
  # so either that object or the bare vector is accepted here
  if (!is.null(reference) && length(reference) == 1L &&
      !is.null(attr(reference, "reference")))
    reference <- attr(reference, "reference")
  d <- as.data.frame(x)
  if (!is.null(reference) && length(reference) < nrow(d))
    stop("`reference` must have one value per component (", nrow(d), ").",
         call. = FALSE)
  d$component <- factor(d$component, levels = d$component)
  plot <- ggplot2::ggplot(d, ggplot2::aes(.data$component, .data$eigenvalue,
                                          group = 1)) +
    ggplot2::geom_line(colour = "grey40") +
    ggplot2::geom_point(size = 2, colour = "#2166AC") +
    ggplot2::labs(x = "Component", y = "Eigenvalue", title = title) +
    ggplot2::theme_minimal(base_size = 11)
  if (!is.null(reference)) {
    line <- data.frame(component = d$component,
                       eigenvalue = reference[seq_len(nrow(d))])
    plot <- plot + ggplot2::geom_line(data = line, linetype = "dashed",
                                      colour = "#B2182B")
  }
  plot
}

# Eigenvalues of the correlation matrix of the residuals, descending. Traits
# are checked before this point, but a trait can still be fitted exactly by the
# fixed effects, which leaves nothing to standardize. The tolerance is relative
# because an exact fit leaves round-off rather than a zero.
.residual_eigenvalues <- function(Y, X, traits = NULL) {
  residuals <- if (ncol(X) == 1L && all(X == 1)) Y else stats::lm.fit(X, Y)$residuals
  spread <- apply(residuals, 2L, stats::sd)
  dead <- spread <= 1e-8 * apply(Y, 2L, stats::sd)
  if (any(dead)) {
    labels <- traits %||% paste("column", seq_len(ncol(Y)))
    stop("`X` fits ", paste(labels[dead], collapse = ", "),
         " exactly, so there is no residual variance to decompose.",
         call. = FALSE)
  }
  eigen(stats::cor(residuals), symmetric = TRUE, only.values = TRUE)$values
}

.scree_frame <- function(values) {
  total <- sum(values)
  out <- data.frame(component = seq_along(values),
                    eigenvalue = values,
                    share = values / total,
                    cumulative = cumsum(values) / total)
  structure(out, class = c("scree", "data.frame"))
}
