# Compare the five covariance structures available in lmm_animal()
#
# Simulates K = 5 traits from a two-factor genetic covariance G and a diagonal
# residual covariance E, then fits the five multivariate models available in
# stanimal:
#
#   1. unstructured G, unstructured E   (the default)
#   2. unstructured G, diagonal E
#   3. factor G,       unstructured E
#   4. factor G,       diagonal E       (the true model here)
#   5. factor G,       shared-factor E
#
# The shared-factor model uses the same latent factors for G and E; with a
# truly diagonal E it should learn factor heritabilities near one. The
# comparison follows examples/trait_comparison.R: relative MSE of G, E,
# heritability and genetic correlations, the minimum bulk ESS of G, sampler
# divergences, and the time each fit took.
#
# Relative MSE scores a model against the values that generated the data, which
# is only available in a simulation. The last part of the script adds the
# question a user can actually ask of real data: which structure predicts a
# record it has not seen? Every model reports a pointwise log_lik, so fit_loo()
# gives each one an expected log predictive density, and loo_compare() ranks
# them with the standard error of the difference between two of them. Run from
# the project root:
#   Rscript examples/structure_comparison.R

if (!requireNamespace("cmdstanr", quietly = TRUE))
  stop("Install cmdstanr and CmdStan first: https://mc-stan.org/cmdstanr/")
if (!requireNamespace("loo", quietly = TRUE))
  stop("Install the loo package first: install.packages(\"loo\").")

library(stanimal)

## Simulation ----------------------------------------------------------------

set.seed(1234)

n_traits <- 5L
n_factors <- 2L
traits <- paste0("t", seq_len(n_traits))

# A sparse two-factor genetic covariance: the first factor loads on traits 1-3,
# the second on traits 3-5. The diagonal Psi keeps G positive definite.
Lambda <- matrix(c(
  0.90,  0.00,
  0.85,  0.00,
  0.70,  0.30,
  0.00,  0.80,
  0.00, -0.70
), n_traits, n_factors, byrow = TRUE)
Psi <- c(0.35, 0.40, 0.30, 0.45, 0.50)
G <- Lambda %*% t(Lambda) + diag(Psi)
E <- diag(c(0.80, 0.90, 0.75, 0.85, 0.95))
colnames(G) <- rownames(G) <- traits
colnames(E) <- rownames(E) <- traits

# Five generations of 100 individuals, with sex-consistent parents.
n_gen <- 5; n_breed <- 100
ped <- data.frame(id = seq_len(n_gen * n_breed),
                  sire = NA_integer_, dam = NA_integer_,
                  sex = NA_character_, stringsAsFactors = FALSE)
ped$sex[seq_len(n_breed)] <- sample(rep(c("M", "F"), each = n_breed / 2))
for (g in 2:n_gen) {
  kids <- ((g - 1) * n_breed + 1):(g * n_breed)
  pool <- ped[ped$id < min(kids), ]
  ped$sire[kids] <- sample(pool$id[pool$sex == "M"], n_breed, TRUE)
  ped$dam[kids]  <- sample(pool$id[pool$sex == "F"], n_breed, TRUE)
  ped$sex[kids]  <- sample(c("M", "F"), n_breed, TRUE)
}
A <- makeA(ped)

n <- nrow(A)
X <- cbind(Intercept = 1, sex = rbinom(n, 1, 0.5), age = round(rnorm(n, 3, 0.5), 2))
B <- rbind(rep(10, n_traits), rep(1.5, n_traits), rep(0.4, n_traits))

Y <- X %*% B +
  t(chol(A)) %*% matrix(rnorm(n * n_traits), n) %*% chol(G) +
  matrix(rnorm(n * n_traits), n) %*% chol(E)
colnames(Y) <- traits

## Comparison helpers --------------------------------------------------------

corr_of <- function(S) {
  if (nrow(S) < 2L) return(numeric(0))
  cov2cor(S)[upper.tri(S)]
}

h2_true <- diag(G) / diag(G + E)
rG_true <- corr_of(G)

# Relative MSE makes errors comparable across quantities with different scales.
relative_mse <- function(estimate, truth)
  mean((estimate - truth)^2) / mean(truth^2)

# One row per model: relative MSE of G, E, the heritabilities and the genetic
# correlations, held-out predictive density, the minimum ESS of G, ESS per
# retained draw, and time.
compare <- function(model, G_hat, E_hat, elapsed = NA_real_,
                    min_ess = NA_real_, n_draws = NA_real_,
                    divergences = NA_real_, elpd = NA_real_,
                    elpd_se = NA_real_, p_loo = NA_real_,
                    pareto_k_high = NA_real_) {
  ok <- !is.null(G_hat) && !is.null(E_hat) &&
    !anyNA(G_hat) && !anyNA(E_hat)
  if (!ok)
    return(data.frame(model = model, G_relative_mse = NA_real_,
                      E_relative_mse = NA_real_,
                      h2_relative_mse = NA_real_,
                      corrG_relative_mse = NA_real_,
                      elpd = elpd, se_elpd = elpd_se, p_loo = p_loo,
                      pareto_k_high = pareto_k_high, min_ess = min_ess,
                      ess_ratio = min_ess / n_draws,
                      divergences = divergences, time_s = elapsed))
  data.frame(model = model,
             G_relative_mse = relative_mse(G_hat, G),
             E_relative_mse = relative_mse(E_hat, E),
             h2_relative_mse = relative_mse(
               diag(G_hat) / diag(G_hat + E_hat), h2_true),
             corrG_relative_mse = relative_mse(corr_of(G_hat), rG_true),
             elpd = elpd, se_elpd = elpd_se, p_loo = p_loo,
             pareto_k_high = pareto_k_high,
             min_ess = min_ess, ess_ratio = min_ess / n_draws,
             divergences = divergences, time_s = elapsed)
}

# Smallest bulk ESS across the unique elements of G.
ess_min <- function(draws) {
  draws <- as.matrix(draws)
  if (!nrow(draws) || !ncol(draws)) return(NA_real_)
  values <- vapply(seq_len(ncol(draws)),
                   function(j) posterior::ess_bulk(draws[, j]), numeric(1))
  if (all(is.na(values))) return(NA_real_)
  min(values, na.rm = TRUE)
}

# Fit one model, then compute the same summaries for every row of the table.
# The LOO object is kept whole so that loo_compare() can use it below.
fit_model <- function(label, ...) {
  elapsed <- system.time(
    fit <- lmm_animal(Y, X, A, iter = 2000, warmup = 1000,
                      chains = 4, cores = 4, refresh = 0, ...)
  )[["elapsed"]]

  g_draws <- fit_draws(fit)$G
  g_pairs <- which(lower.tri(diag(n_traits), diag = TRUE), arr.ind = TRUE)
  g_matrix <- vapply(seq_len(nrow(g_pairs)), function(i)
    g_draws[g_pairs[i, 1], g_pairs[i, 2], ], numeric(dim(g_draws)[3]))

  diagnostics <- fit$diagnostic_summary(quiet = TRUE)
  # the Pareto-k warning would repeat once per model; k above the threshold is a
  # column of the table below, which is more use than five copies of it
  loo_fit <- withCallingHandlers(
    fit_loo(fit),
    warning = function(w) {
      if (grepl("Pareto k diagnostic", conditionMessage(w)))
        invokeRestart("muffleWarning")
    })
  estimates <- loo_fit$estimates
  list(
    row = compare(label,
                  draws_mean(fit, "G"), draws_mean(fit, "E"),
                  elapsed = elapsed,
                  min_ess = ess_min(g_matrix),
                  n_draws = nrow(g_matrix),
                  divergences = sum(diagnostics$num_divergent),
                  elpd = unname(estimates["elpd_loo", "Estimate"]),
                  elpd_se = unname(estimates["elpd_loo", "SE"]),
                  p_loo = unname(estimates["p_loo", "Estimate"]),
                  # k above 0.7 means the importance weights are unreliable for
                  # that record, which is worth seeing next to the total
                  pareto_k_high = sum(loo_fit$diagnostics$pareto_k > 0.7)),
    loo = loo_fit)
}

## The five models -----------------------------------------------------------

fits <- list()

fits[["G us, E us"]] <- fit_model(
  "G us, E us",
  residual_structure = "unstructured",
  genetic_structure = "unstructured"
)

fits[["G us, E diag"]] <- fit_model(
  "G us, E diag",
  residual_structure = "diagonal",
  genetic_structure = "unstructured"
)

fits[["G factor, E us"]] <- fit_model(
  "G factor, E us",
  residual_structure = "unstructured",
  genetic_structure = "factor", n_factors = n_factors,
  hs_global_scale = 2, hs_slab_scale = 2
)

fits[["G factor, E diag"]] <- fit_model(
  "G factor, E diag",
  residual_structure = "diagonal",
  genetic_structure = "factor", n_factors = n_factors,
  hs_global_scale = 2, hs_slab_scale = 2
)

fits[["G factor, E shared"]] <- fit_model(
  "G factor, E shared",
  residual_structure = "shared_factor",
  genetic_structure = "factor", n_factors = n_factors,
  hs_global_scale = 2, hs_slab_scale = 2
)

## The comparison ------------------------------------------------------------

results <- rbind(compare("truth", G, E),
                 do.call(rbind, lapply(fits, `[[`, "row")))
rownames(results) <- results$model
out <- round(results[, -1], 4)
out[, "time_s"] <- round(results$time_s, 1)
out[, "divergences"] <- as.integer(results$divergences)
out[, "pareto_k_high"] <- as.integer(results$pareto_k_high)
cat("\n", n_traits, "traits,", n, "individuals,", n_factors,
    "true genetic factors \n", sep = "")
print(out)

# The ranking, with the standard error of each difference from the best model.
# A difference smaller than about twice its standard error is not worth reading
# as an ordering. loo_compare() takes the list of LOO objects and labels each row
# with the name it had there; recent versions return a data frame with a `model`
# column, earlier ones a matrix whose rownames are the models.
comparison <- as.data.frame(loo::loo_compare(lapply(fits, `[[`, "loo")))
if (!"model" %in% names(comparison)) {
  comparison$model <- rownames(comparison)
  rownames(comparison) <- NULL
}
shown <- intersect(c("model", "elpd_diff", "se_diff", "elpd_loo", "p_loo"),
                   names(comparison))
comparison <- comparison[, shown, drop = FALSE]
numeric_columns <- vapply(comparison, is.numeric, logical(1))
comparison[numeric_columns] <- lapply(comparison[numeric_columns], round, 1)
cat("\nHeld-out predictive density, best first (higher ELPD is better):\n")
print(comparison, row.names = FALSE)

# Relative MSE divides each elementwise mean squared error by the mean squared
# truth on the same scale. The correlation comparison uses the upper triangle
# because its diagonal is always 1. min_ess is the smallest bulk ESS across the
# unique elements of G, and ess_ratio divides it by the 4000 retained draws
# (4 chains x 1000), so it measures mixing rather than chain length. The truth
# row has no draws, no time and no ELPD.
#
# `elpd` is the expected log pointwise predictive density from PSIS-LOO, and
# `p_loo` its effective number of parameters, so a structure that buys its fit
# with many parameters shows it here. `pareto_k_high` counts the records whose
# importance weights the Pareto fit cannot be trusted for (k > 0.7): the total
# ELPD is only as good as those. In the loo_compare() table below, every model
# is compared with the best one, so `elpd_diff` is its loss relative to that
# model and `se_diff` is the uncertainty of that difference; a difference well
# inside a couple of standard errors is not an ordering.
#
# What the ranking separates, and what it does not, is the interesting part. The
# residual structure it does separate: both unstructured-E models lose 6.5 to 7.7
# ELPD against standard errors near 2.4, so the data prefer a diagonal E, which is
# what generated them. The genetic structure it does not: within either residual
# choice the factor and unstructured G are 1.3 ELPD apart at most, against
# standard errors of 1.0 to 1.4. If the factor structure is to be preferred here,
# the argument is `p_loo` (32 effective parameters against 42) rather than the
# held-out density, and that argument gets stronger as the number of traits grows,
# which is where the factor model stops being a rounding error.
#
# Two things keep the ELPD comparison honest. Every fit is to the same response,
# so the constant that the standardization contributes to each term cancels.
# And the terms are those of the model block, one per record of the *rotated*
# response: their sum is the exact joint log density, but with a non-identity A
# a term is not one individual left out. That is why this compares structures
# and not individuals -- for the individual-level question, use
# posterior_predictive().
#
# Timings are a single run each, in seconds, and cover the lmm_animal() call
# only. The first stanimal fit may include compilation when the model cache is
# cold; later fits reuse the compiled executable. The unstructured models use
# the same LKJ priors, and the factor models use the same horseshoe and column
# priors, so differences reflect the covariance structure rather than different
# sampling settings. With only five traits the factor model has
# little parameter savings, and its hierarchical shrinkage can cost a little
# accuracy and time; its advantage is expected to grow with the number of
# traits, where an unstructured G becomes expensive.
