# Completed datasets from a fit with incomplete records.
#
# Both engines integrate the missing values out, so neither reports a completed
# response as a by-product. What they do give is enough to draw one:
#
#   - the animal models sample the gaps as parameters (`Y_mis`, on the
#     standardized response scale), so each draw already is an imputation, and
#     undoing the standardization recovers it on the data scale;
#   - the phenotypic models have no gap parameters, because a row's missing
#     values can be integrated out analytically. There the gap is drawn from its
#     conditional distribution given the row's observed traits, which is the
#     same thing: the posterior of a missing value under that model.
#
# Either way the result is *multiple* imputation -- one completed dataset per
# draw -- not a single filled-in matrix, which would understate the uncertainty
# in exactly the cells nobody measured.

# Reference values for the standardization the model block applies: the mean and
# standard deviation of the observed entries of each trait, the latter with the
# n - 1 divisor, matching transformed data in every Stan file.
.standardization <- function(Y) {
  miss <- is.na(Y)
  lapply(seq_len(ncol(Y)), function(j) {
    observed <- Y[!miss[, j], j]
    list(mean = mean(observed), sd = stats::sd(observed))
  })
}

# Which posterior draws to use, spread over the kept draws so that both chains
# are represented.
.draw_index <- function(kept, ndraws) {
  ndraws <- as.integer(ndraws)
  if (ndraws < 1L)
    stop("`ndraws` must be at least 1.", call. = FALSE)
  if (ndraws >= kept) seq_len(kept) else unique(round(seq(1, kept, length.out = ndraws)))
}

#' Completed datasets from a fit with missing values
#'
#' Draws the response with its gaps filled in, one completed dataset per
#' posterior draw. The observed entries are the data and are the same in every
#' draw; only the cells that were not measured vary, which is how the uncertainty
#' about them shows up.
#'
#' @param fit A fit from \code{\link{lmm_animal}} or \code{\link{lm_phen_cov}}
#'   whose response had missing values.
#' @param ndraws Number of completed datasets to draw, taken evenly across the
#'   kept draws.
#' @param seed Optional seed for the draws that are simulated rather than read
#'   from the sampler. It only matters for a phenotypic fit; an animal fit's gaps
#'   are already parameters.
#' @return An object of class \code{stanimal_imputation}: a list with
#'   \code{Y} (individuals by traits by draws, with draws last, as in
#'   \code{\link{fit_draws}}), the \code{observed} response it came from, the
#'   \code{gaps} as a data frame of individual and trait positions, and the trait
#'   names. \code{print()} reports the imputed value of each gap with its 95\%
#'   interval.
#' @details An animal fit's gaps are parameters, so each draw is already a draw
#'   from the posterior of the missing values and is reported as it is. A
#'   phenotypic fit has no such parameters -- the missing values were integrated
#'   out -- so each gap is drawn from the normal conditional distribution given
#'   that row's observed traits, using \code{P} and \code{beta} at the same draw.
#'   The two give the same kind of object, which is why one function covers both.
#'
#'   Use these for imputation-style follow-up analyses, or to see what the model
#'   believes about the cells it never saw. Do not collapse them to a single
#'   filled-in matrix and treat it as data: the spread across draws is the point.
#' @seealso \code{\link{posterior_predictive}}, \code{\link{fit_loo}}
#' @examples
#' \dontrun{
#' Y[c(2, 9), 1] <- NA
#' fit <- lmm_animal(Y, X, A)
#' imp <- imputed_data(fit, ndraws = 200)
#' imp                                   # each gap with its interval
#' dim(imp$Y)                            # individuals x traits x draws
#' }
#' @export
imputed_data <- function(fit, ndraws = 100, seed = NULL) {
  data <- .stanimal_data(fit)
  Y <- data$Y
  miss <- is.na(Y)
  if (!any(miss))
    stop("the response this fit was given has no missing values, so there is ",
         "nothing to impute.", call. = FALSE)
  n <- nrow(Y)
  k <- ncol(Y)
  traits <- .draw_traits(fit, NULL, k)
  draws <- .draws_matrix(fit)
  take <- .draw_index(nrow(draws), ndraws)

  # Positions of the gaps in column-major order, which is the order the gap
  # parameters were declared in, and the order R assigns them below.
  gaps <- data.frame(individual = row(miss)[miss], trait = col(miss)[miss])
  gaps$individual_label <- rownames(Y)[gaps$individual] %||%
    as.character(gaps$individual)
  gaps$trait_label <- traits[gaps$trait]
  standard <- .standardization(Y)

  imputed <- if (data$family == "animal") {
    .impute_animal(draws, take, gaps, standard)
  } else {
    if (!is.null(seed)) set.seed(seed)
    .impute_phenotypic(draws, take, Y, data$X, gaps, k, traits, fit)
  }

  completed <- array(NA_real_, c(n, k, length(take)),
                     dimnames = list(rownames(Y), traits, NULL))
  for (s in seq_along(take)) {
    value <- Y
    value[miss] <- imputed[s, ]
    completed[, , s] <- value
  }
  structure(list(Y = completed, observed = Y, gaps = gaps, traits = traits),
            class = "stanimal_imputation")
}

# The gaps the animal model sampled, undone back to the response scale.
.impute_animal <- function(draws, take, gaps, standard) {
  columns <- grep("^Y_mis", colnames(draws))
  if (!length(columns))
    stop("this animal fit does not report `Y_mis`, so the gaps cannot be ",
         "recovered from its draws.", call. = FALSE)
  if (length(columns) != nrow(gaps))
    stop("this fit reports ", length(columns), " gap parameters but the ",
         "response has ", nrow(gaps), " gaps.", call. = FALSE)
  values <- draws[take, columns, drop = FALSE]
  mean_scale <- vapply(standard[gaps$trait], `[[`, numeric(1), "mean")
  sd_scale <- vapply(standard[gaps$trait], `[[`, numeric(1), "sd")
  sweep(sweep(values, 2L, sd_scale, "*"), 2L, mean_scale, "+")
}

# The gaps the phenotypic model integrated out, drawn from their conditional
# distribution given the observed traits of the same row.
.impute_phenotypic <- function(draws, take, Y, X, gaps, k, traits, fit) {
  covariance <- .matrix_draws(draws, "P", k, traits)
  beta <- .beta_draws(draws, k, traits, NULL, attr(fit, "x_names"))
  if (is.null(covariance) || is.null(beta))
    stop("this phenotypic fit does not report `P` and `beta`, so the gaps ",
         "cannot be drawn.", call. = FALSE)
  miss <- is.na(Y)
  key <- apply(!miss, 1L, function(o) paste(as.integer(o), collapse = ""))
  out <- matrix(NA_real_, nrow = length(take), ncol = nrow(gaps))
  for (s in seq_along(take)) {
    j <- take[s]
    Sig <- matrix(covariance[, , j], k, k)
    mu <- X %*% t(matrix(beta[, , j], nrow = k))
    for (pattern in unique(key)) {
      rows <- which(key == pattern)
      observed <- which(!miss[rows[1L], ])
      missing <- setdiff(seq_len(k), observed)
      if (!length(missing)) next
      # Sigma over the observed traits, jittered for the draws where a factor
      # structure leaves it numerically singular
      Soo <- Sig[observed, observed, drop = FALSE] +
        .jitter(Sig[observed, observed, drop = FALSE])
      Smo <- Sig[missing, observed, drop = FALSE]
      B <- Smo %*% solve(Soo)
      cond_cov <- Sig[missing, missing, drop = FALSE] - B %*% t(Smo)
      cond_cov <- (cond_cov + t(cond_cov)) / 2 + .jitter(cond_cov)
      L <- chol(cond_cov)
      resid <- Y[rows, observed, drop = FALSE] - mu[rows, observed, drop = FALSE]
      cond_mean <- mu[rows, missing, drop = FALSE] + resid %*% t(B)
      noise <- matrix(stats::rnorm(length(rows) * length(missing)),
                      length(rows), length(missing)) %*% L
      filled <- cond_mean + noise
      # place the drawn values in the gap columns of this draw
      wanted <- which(gaps$trait %in% missing)
      for (r in seq_along(rows)) {
        row_gaps <- which(gaps$individual == rows[r])
        if (!length(row_gaps)) next
        out[s, row_gaps] <- filled[r, match(gaps$trait[row_gaps], missing)]
      }
    }
  }
  out
}

#' @rdname imputed_data
#' @param x An object from \code{imputed_data()}.
#' @param ... Unused; for method compatibility.
#' @export
print.stanimal_imputation <- function(x, ...) {
  gaps <- x$gaps
  draws <- x$Y[, , , drop = FALSE]
  values <- t(vapply(seq_len(nrow(gaps)), function(i) {
    v <- draws[gaps$individual[i], gaps$trait[i], ]
    c(mean = mean(v),
      q2.5 = stats::quantile(v, 0.025, names = FALSE),
      q97.5 = stats::quantile(v, 0.975, names = FALSE))
  }, numeric(3)))
  out <- data.frame(individual = gaps$individual_label,
                    trait = gaps$trait_label,
                    round(values, 3), check.names = FALSE)
  # the gaps are stored in the column-major order the parameters were declared
  # in, which is not the order a person reads them in
  out <- out[order(gaps$individual, gaps$trait), , drop = FALSE]
  rownames(out) <- NULL
  cat("<stanimal_imputation>", dim(draws)[3L], "completed datasets of",
      dim(draws)[1L], "individuals and", dim(draws)[2L], "trait(s),",
      nrow(gaps), "imputed value(s)\n")
  print(out, row.names = FALSE)
  invisible(x)
}
