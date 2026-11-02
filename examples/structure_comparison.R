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
# comparison follows examples/trait_comparison.R: elementwise MSE of G, E,
# heritability and genetic correlations, the minimum bulk ESS of G, sampler
# divergences, and the time each fit took. Run from the project root:
#   Rscript examples/structure_comparison.R

if (!requireNamespace("cmdstanr", quietly = TRUE))
  stop("Install cmdstanr and CmdStan first: https://mc-stan.org/cmdstanr/")

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

# One row per model: MSE of G, E, the heritabilities and the genetic
# correlations, the minimum ESS of G, ESS per retained draw, and time.
compare <- function(model, G_hat, E_hat, elapsed = NA_real_,
                    min_ess = NA_real_, n_draws = NA_real_,
                    divergences = NA_real_) {
  ok <- !is.null(G_hat) && !is.null(E_hat) &&
    !anyNA(G_hat) && !anyNA(E_hat)
  if (!ok)
    return(data.frame(model = model, G_mse = NA_real_, E_mse = NA_real_,
                      h2_mse = NA_real_, corrG_mse = NA_real_,
                      min_ess = min_ess, ess_ratio = min_ess / n_draws,
                      divergences = divergences, time_s = elapsed))
  data.frame(model = model,
             G_mse = mean((G_hat - G)^2),
             E_mse = mean((E_hat - E)^2),
             h2_mse = mean((diag(G_hat) / diag(G_hat + E_hat) - h2_true)^2),
             corrG_mse = mean((corr_of(G_hat) - rG_true)^2),
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
  compare(label,
          draws_mean(fit, "G"), draws_mean(fit, "E"),
          elapsed = elapsed,
          min_ess = ess_min(g_matrix),
          n_draws = nrow(g_matrix),
          divergences = sum(diagnostics$num_divergent))
}

## The five models -----------------------------------------------------------

results <- list(compare("truth", G, E))

results[["G us, E us"]] <- fit_model(
  "G us, E us",
  residual_structure = "unstructured",
  genetic_structure = "unstructured"
)

results[["G us, E diag"]] <- fit_model(
  "G us, E diag",
  residual_structure = "diagonal",
  genetic_structure = "unstructured"
)

results[["G factor, E us"]] <- fit_model(
  "G factor, E us",
  residual_structure = "unstructured",
  genetic_structure = "factor", n_factors = n_factors,
  hs_global_scale = 2, hs_slab_scale = 2
)

results[["G factor, E diag"]] <- fit_model(
  "G factor, E diag",
  residual_structure = "diagonal",
  genetic_structure = "factor", n_factors = n_factors,
  hs_global_scale = 2, hs_slab_scale = 2
)

results[["G factor, E shared"]] <- fit_model(
  "G factor, E shared",
  residual_structure = "shared_factor",
  genetic_structure = "factor", n_factors = n_factors,
  hs_global_scale = 2, hs_slab_scale = 2
)

## The comparison ------------------------------------------------------------

results <- do.call(rbind, results)
rownames(results) <- results$model
out <- round(results[, -1], 4)
out[, "time_s"] <- round(results$time_s, 1)
out[, "divergences"] <- as.integer(results$divergences)
cat("\n", n_traits, "traits,", n, "individuals,", n_factors,
    "true genetic factors \n", sep = "")
print(out)

# G_mse and E_mse are elementwise mean squared errors against the simulated
# matrices. h2_mse and corrG_mse compare the trait-level summaries, with corrG
# using the upper triangle because its diagonal is always 1. min_ess is the
# smallest bulk ESS across the unique elements of G, and ess_ratio divides it
# by the 4000 retained draws (4 chains x 1000), so it measures mixing rather
# than chain length. The truth row has no draws and no time.
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
