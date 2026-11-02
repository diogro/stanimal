# Sparse factor model for a modular genetic covariance
#
# Simulates twelve traits from a nine-factor hierarchy: a general size factor,
# two broad blocks, four asymmetric sub-blocks, and two nested sub-sub-blocks.
# It fits genetic_structure = "factor" and compares the posterior G, genetic
# correlations, h2, uniqueness and factor magnitudes with the simulation. Run
# from the project root:
#   Rscript examples/factor_model.R

if(!requireNamespace("cmdstanr", quietly = TRUE))
  stop("Install cmdstanr and CmdStan first: https://mc-stan.org/cmdstanr/")
library(stanimal)

set.seed(1)
n_traits <- 12L
n_factors <- 9L
traits <- sprintf("trait%d", seq_len(n_traits))

# A hierarchical loading matrix. The general size factor loads on every trait;
# the two blocks split the traits 1-6 and 7-12; the sub-blocks are deliberately
# unequal (2, 4, 3 and 3 traits); and only two of the sub-blocks contain a
# sub-sub-block, so the matrix is heterogeneous rather than a symmetric grid.
Lambda <- matrix(0, n_traits, n_factors)
colnames(Lambda) <- c("size", "block1", "block2",
                      "sub1a", "sub1b", "sub2a", "sub2b",
                      "subsub1", "subsub2")
Lambda[, 1] <- c(0.62, 0.58, 0.55, 0.50, 0.48, 0.52,
                 0.60, 0.56, 0.53, 0.49, 0.57, 0.51)
Lambda[1:6, 2] <- c(0.55, 0.50, 0.45, 0.48, 0.52, 0.40)
Lambda[7:12, 3] <- c(0.50, 0.58, 0.45, 0.55, 0.48, 0.52)
Lambda[1:2, 4] <- c(0.55, 0.50)       # sub-block 1a: traits 1-2
Lambda[3:6, 5] <- c(0.50, 0.45, 0.55, 0.40)  # sub-block 1b: traits 3-6
Lambda[7:9, 6] <- c(0.45, 0.50, 0.40) # sub-block 2a: traits 7-9
Lambda[10:12, 7] <- c(0.55, 0.50, 0.45) # sub-block 2b: traits 10-12
Lambda[3:4, 8] <- c(0.45, 0.40)       # nested in sub-block 1b
Lambda[10:11, 9] <- c(0.45, 0.40)     # nested in sub-block 2b

Psi <- c(0.35, 0.40, 0.30, 0.45, 0.50, 0.35,
         0.40, 0.30, 0.45, 0.35, 0.40, 0.45)
G <- Lambda %*% t(Lambda) + diag(Psi)
E <- diag(c(0.80, 0.90, 0.75, 0.85, 0.95, 0.80,
            0.90, 0.75, 0.85, 0.80, 0.90, 0.75))
rownames(G) <- colnames(G) <- traits
rownames(E) <- colnames(E) <- traits

# Five generations of 120 individuals, with sex-consistent parents.
n_gen <- 5L; n_breed <- 120L
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
X <- cbind(Intercept = 1, sex = rbinom(n, 1, 0.5),
           age = round(rnorm(n, 3, 0.5), 2))
B <- rbind(rep(10, n_traits),
           seq(-1, 1.5, length.out = n_traits),
           seq(-0.4, 0.7, length.out = n_traits))
u <- t(chol(A)) %*% matrix(rnorm(n * n_traits), n) %*% chol(G)
e <- matrix(rnorm(n * n_traits), n) %*% chol(E)
Y <- X %*% B + u + e
colnames(Y) <- traits

fit <- lmm_animal(
  Y, X, A,
  genetic_structure = "factor", n_factors = n_factors,
  hs_global_scale = 2, hs_slab_scale = 2,
  # The hierarchy is ordered from broad to narrow; a mild decay favours that
  # order without forcing it.
  column_decay = 0.9,
  iter = 2500, warmup = 1250, chains = 4, cores = 4,
  refresh = 0, show_messages = FALSE
)

# The posterior mean of G is close to the simulated factor covariance. The
# error is reported relative to the size of G, so it stays readable whatever
# scale the simulation uses; an absolute RMSE on a covariance matrix means
# nothing without the variances it should be compared with.
G_hat <- draws_mean(fit, "G")
cat("G relative error:",
    round(sqrt(mean((G_hat - G)^2)) / sqrt(mean(G^2)), 3), "\n")
round(G_hat, 2)
round(G, 2)

# Genetic correlations recover the overlapping module structure.
round(draws_mean(fit, "corrG"), 2)
round(cov2cor(G), 2)

# Heritabilities and communalities are the easiest summaries to compare with
# the simulation. Loadings themselves are only identified up to rotation and
# sign, so their posterior means can be near zero even when G is recovered.
h2_hat <- draws_mean(fit, "h2")
h2_true <- diag(G) / diag(G + E)
cat("h2 RMSE:", round(sqrt(mean((h2_hat - h2_true)^2)), 3), "\n")
round(h2_hat, 2)
round(h2_true, 2)

communality_hat <- draws_mean(fit, "communality")
communality_true <- rowSums(Lambda^2) / diag(G)
round(communality_hat, 2)
round(communality_true, 2)
round(draws_mean(fit, "uniqueness"), 2)
round(Psi / diag(G), 2)

# factor_magnitude is the squared norm of each loading column on the
# correlation scale. Compare it with the simulated column norms.
fit$summary("factor_magnitude")[, c("variable", "mean", "sd", "rhat", "ess_bulk")]
round(colSums(Lambda^2), 2)

# Sampler diagnostics. The pointwise log likelihood is also available as
# fit$draws("log_lik") for model comparison with loo.
print(fit$diagnostic_summary())

# Heatmaps of the estimated and simulated genetic correlations.
dir.create("examples/figures", recursive = TRUE, showWarnings = FALSE)
pdf("examples/figures/factor_model_corr.pdf", width = 11, height = 5)
print(plot_corr(draws_mean(fit, "corrG"), "Posterior mean genetic correlation"))
print(plot_corr(cov2cor(G), "Simulated genetic correlation"))
dev.off()
