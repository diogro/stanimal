# Posterior predictive replicates of the response.
#
# `log_lik` answers "how well does this model predict a record it has not seen",
# which is what `fit_loo()` needs but is an abstract number. Replicating the
# response from the posterior answers the question a geneticist actually asks:
# if the fitted model were true, would the data look like the data?
#
# The replicate is drawn from the marginal model the fit describes, in the
# original response scale: the fixed effects plus the genetic and residual
# deviations (animal fits) or the residual deviations alone (phenotypic fits).
# Only quantities the fit reports are used, so factor and unstructured fits are
# treated alike.

# `replicated` is a column of the data frame the plot method builds; R CMD check
# cannot see through ggplot2's non-standard evaluation.
utils::globalVariables("replicated")

# The response, design and relationship matrix the fit was given.
.stanimal_data <- function(fit, what = "the response cannot be replicated") {
  data <- attr(fit, "stanimal_data")
  if (is.null(data))
    stop("this fit does not carry the data it was fitted to, so ", what, ". ",
         "Fits from `lmm_animal()` and `lm_phen_cov()` do; ",
         "`lmm_animal_gibbs()` and other objects do not.", call. = FALSE)
  data
}

#' Replicate the response from the posterior
#'
#' Draws replicate datasets from the fitted model, one per posterior draw, so
#' that the data can be compared with what the model says they should look like.
#' This is the individual-level check that \code{\link{fit_loo}} is not: LOO
#' compares model structures by held-out predictive density, while this shows
#' \emph{which} features of the data the fit reproduces.
#'
#' @param fit A fit from \code{\link{lmm_animal}} or \code{\link{lm_phen_cov}}.
#' @param ndraws Number of posterior draws to replicate from; taken evenly across
#'   the kept draws, so both warmup-free chains are represented.
#' @param seed Optional seed for the simulation, to make the replicates
#'   reproducible. The sampler's seed does not cover this step.
#' @param ... Unused; for method compatibility.
#' @return An object of class \code{stanimal_yrep}: a list with
#'   \code{yrep} (individuals by traits by draws, the same layout as
#'   \code{\link{fit_draws}}, with draws last), the observed \code{y} on the
#'   original response scale, the trait names, and the number of replicates.
#'   \code{print()} reports, for the trait means, standard deviations and
#'   pairwise correlations, the observed value, the 95\% interval of the
#'   replicates, and \code{p}: the share of replicates at or above the observed
#'   value, a posterior predictive tail probability. These probabilities are
#'   not generally uniform under a correct model because the data are used
#'   both to fit and to check it. Values near 0 or 1 can flag a discrepancy;
#'   a middling value does not establish a good fit. \code{plot()} shows the
#'   replicate distributions.
#' @details For an animal fit the replicate is
#'   \eqn{X\beta + u + e} with \eqn{u \sim N(0, A \otimes G)} and
#'   \eqn{e \sim N(0, I \otimes E)}; for a phenotypic fit it is
#'   \eqn{X\beta + e} with \eqn{e \sim N(0, P)}. \code{G}, \code{E} and \code{P}
#'   are the reported draws, so they are already on the response scale.
#'
#'   Replicates are simulated for every trait of every individual, including the
#'   entries that were not observed, so \code{y} is the response as it was passed
#'   in and has \code{NA} where the data are incomplete. Printed and plotted
#'   summaries compare each replicate on the same observed entries as the
#'   corresponding statistic in \code{y}. Pairs with fewer than two jointly
#'   observed values or constant observed values are omitted with a warning.
#' @seealso \code{\link{fit_loo}}
#' @examples
#' \dontrun{
#' pp <- posterior_predictive(fit)
#' pp                      # observed against replicated summaries
#' plot(pp)                # and the distributions behind them
#' }
#' @export
posterior_predictive <- function(fit, ndraws = 100, seed = NULL, ...) {
  data <- .stanimal_data(fit)
  if (!is.null(seed))
    set.seed(seed)
  Y <- data$Y
  X <- data$X
  n <- nrow(Y)
  k <- ncol(Y)
  draws <- .draws_matrix(fit)
  traits <- .draw_traits(fit, NULL, k)

  beta <- .beta_draws(draws, k, traits, NULL, attr(fit, "x_names"))
  if (is.null(beta))
    stop("this fit does not report `beta`, so the response cannot be ",
         "replicated.", call. = FALSE)

  if (data$family == "phenotypic") {
    # the residual covariance is the whole story: y* = X beta + e, e ~ N(0, P)
    covariance <- .matrix_draws(draws, "P", k, traits)
    if (is.null(covariance))
      stop("this fit does not report `P`, so the response cannot be replicated.",
           call. = FALSE)
    G <- NULL
  } else {
    covariance <- .as_matrix_draws(.matrix_draws(draws, "E", k, traits),
                                   .vector_draws(draws, "sigma_E", k, traits), k)
    G <- .as_matrix_draws(.matrix_draws(draws, "G", k, traits),
                          .vector_draws(draws, "sigma_G", k, traits), k)
    if (is.null(covariance) || is.null(G))
      stop("this fit does not report the variances, so the response cannot be ",
           "replicated.", call. = FALSE)
    if (is.null(data$A))
      stop("this fit does not carry the relationship matrix, so the genetic ",
           "deviations cannot be simulated.", call. = FALSE)
  }

  ndraws <- as.integer(ndraws)
  if (ndraws < 1L)
    stop("`ndraws` must be at least 1.", call. = FALSE)
  kept <- dim(beta)[3L]
  take <- if (ndraws >= kept) seq_len(kept) else
    unique(round(seq(1, kept, length.out = ndraws)))

  L_A <- if (data$family == "animal") t(chol(data$A + .jitter(data$A))) else NULL
  # individuals, then traits, then draws: the layout fit_draws() uses, so the
  # draws are last in every array the package returns
  yrep <- array(NA_real_, c(n, k, length(take)),
                dimnames = list(rownames(Y), traits, NULL))
  for (s in seq_along(take)) {
    j <- take[s]
    # matrix(): a one-trait fit would otherwise drop the slice to a vector
    value <- X %*% t(matrix(beta[, , j], nrow = k)) +
      .rmvn_chol(n, k, matrix(covariance[, , j], k, k))
    if (data$family == "animal")
      value <- value + L_A %*% .rmvn_chol(n, k, matrix(G[, , j], k, k))
    yrep[, , s] <- value
  }
  structure(list(yrep = yrep, y = Y, traits = traits, ndraws = length(take)),
            class = "stanimal_yrep")
}

# A positive-definite matrix that chol() accepts, for the rare draw whose
# covariance is numerically singular at the fourth decimal place or so.
.jitter <- function(Sigma) {
  scale <- mean(diag(Sigma))
  diag(1e-10 * max(scale, 1e-8), nrow(Sigma))
}

# Centred normal draws with covariance Sigma, as an individuals by traits matrix.
.rmvn_chol <- function(n, k, Sigma) {
  matrix(stats::rnorm(n * k), n, k) %*% chol(Sigma + .jitter(Sigma))
}

# A covariance reported as a matrix when the fit has one, and rebuilt from the
# standard deviations when it only reports those (the univariate models).
.as_matrix_draws <- function(matrices, standard_deviations, k) {  if (!is.null(matrices))
    return(matrices)
  if (is.null(standard_deviations))
    return(NULL)
  draws <- vapply(seq_len(ncol(standard_deviations)),
                  function(s) diag(standard_deviations[, s]^2, nrow = k),
                  if (k == 1L) numeric(1) else matrix(0, k, k))
  array(draws, dim = c(k, k, ncol(standard_deviations)))
}

# The posterior predictive tail probability: the share of replicates at or above
# the observation.
.predictive_tail <- function(replicated, observed) mean(replicated >= observed)

# Observed and replicated values of the summaries that describe a covariance
# structure: trait means, trait standard deviations and pairwise correlations.
#
# Grouped by kind rather than by trait, so all the means sit together and can be
# compared down the column, then all the standard deviations, then all the
# correlations. Both the printed table and the facets of the plot follow this
# order, which is why it is decided here and not in either of them.
.yrep_statistics <- function(x) {
  y <- x$y
  yrep <- x$yrep
  traits <- x$traits
  k <- length(traits)

  # draws are the last dimension, so each statistic is taken over margin 3
  across_draws <- function(f, j) {
    observed <- !is.na(y[, j])
    apply(yrep[observed, j, , drop = FALSE], 3L, f)
  }
  per_trait <- function(prefix, f) {
    do.call(rbind, lapply(seq_len(k), function(j) data.frame(
      statistic = paste0(prefix, traits[j]),
      observed = f(y[, j], na.rm = TRUE),
      replicated = across_draws(f, j))))
  }
  means <- per_trait("mean: ", mean)
  sds <- per_trait("sd: ", stats::sd)

  cors <- list()
  omitted <- character()
  if (k > 1L)
    for (a in seq_len(k - 1L)) for (b in (a + 1L):k) {
      observed <- !is.na(y[, a]) & !is.na(y[, b])
      pair <- paste0(traits[a], "-", traits[b])
      if (sum(observed) < 2L || stats::sd(y[observed, a]) == 0 ||
          stats::sd(y[observed, b]) == 0) {
        omitted <- c(omitted, pair)
        next
      }
      r <- stats::cor(y[observed, a], y[observed, b])
      replicated <- vapply(seq_len(dim(yrep)[3L]), function(s) {
        a_rep <- yrep[observed, a, s]
        b_rep <- yrep[observed, b, s]
        if (any(!is.finite(a_rep)) || any(!is.finite(b_rep)) ||
            stats::sd(a_rep) == 0 || stats::sd(b_rep) == 0)
          stop("replicated correlation is undefined for ", pair, ".",
               call. = FALSE)
        stats::cor(a_rep, b_rep)
      }, numeric(1))
      cors[[length(cors) + 1L]] <- data.frame(
        statistic = paste0("cor: ", pair), observed = r, replicated = replicated)
    }
  if (length(omitted))
    warning("observed correlation is undefined for: ",
            paste(omitted, collapse = ", "), "; omitted from the summaries.",
            call. = FALSE)
  do.call(rbind, c(list(means, sds), cors))
}

#' @rdname posterior_predictive
#' @param x An object returned by \code{posterior_predictive()}.
#' @export
print.stanimal_yrep <- function(x, ...) {
  cat("<stanimal_yrep>", x$ndraws, "replicates of", nrow(x$y), "individuals and",
      length(x$traits), "trait(s)\n")
  s <- .yrep_statistics(x)
  labels <- unique(s$statistic)
  replicated <- lapply(labels, function(g) s$replicated[s$statistic == g])
  observed <- vapply(labels, function(g) s$observed[match(g, s$statistic)],
                     numeric(1))
  out <- data.frame(
    statistic = labels,
    observed = observed,
    q2.5 = vapply(replicated, stats::quantile, numeric(1), probs = 0.025),
    q97.5 = vapply(replicated, stats::quantile, numeric(1), probs = 0.975),
    check.names = FALSE
  )
  names(out)[3:4] <- c("q2.5", "q97.5")
  out$p <- vapply(seq_along(labels), function(i)
    .predictive_tail(replicated[[i]], observed[i]), numeric(1))
  print(out, digits = 3, row.names = FALSE)
  cat("\n`observed` is the data, `q2.5` and `q97.5` are the 95% interval of the\n",
      "replicates, and `p` is the share of replicates at or above the observed\n",
      "value: a posterior predictive tail probability, not a calibrated\n",
      "uniform p-value. Values near 0 or 1 can flag a discrepancy; a middling\n",
      "value does not establish a good fit.\n", sep = "")
  invisible(x)
}

#' @rdname posterior_predictive
#' @param ... Unused; for method compatibility.
#' @export
plot.stanimal_yrep <- function(x, ...) {
  s <- .yrep_statistics(x)
  # the rows of the table are already grouped by kind, so keeping that order as
  # the factor levels puts all the means in the first row of facets, then the
  # standard deviations, then the correlations
  s$statistic <- factor(s$statistic, levels = unique(s$statistic))
  observed <- s[!duplicated(s$statistic), c("statistic", "observed")]
  ggplot2::ggplot(s, ggplot2::aes(x = replicated)) +
    ggplot2::geom_histogram(bins = 30, colour = NA, fill = "grey70") +
    ggplot2::geom_vline(data = observed,
                        ggplot2::aes(xintercept = observed),
                        colour = "firebrick", linewidth = 0.7) +
    ggplot2::facet_wrap(~statistic, scales = "free", ncol = 3) +
    ggplot2::labs(x = "replicated", y = "count of replicates",
                  title = "Posterior predictive check",
                  subtitle = "red line: the observed data") +
    ggplot2::theme_minimal()
}
