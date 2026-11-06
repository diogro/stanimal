# Compare the same phenotypic covariance model with Laplace and sampling.
# Requires stanimal, cmdstanr and a working CmdStan installation.
library(stanimal)

set.seed(2026)
n <- 200
X <- cbind(intercept = 1, predictor = rnorm(n))
P <- matrix(c(1, 0.4, 0.4, 1.5), 2)
Y <- X %*% rbind(c(0.5, -0.5), c(0.3, -0.2)) +
  matrix(rnorm(n * 2), n) %*% chol(P)
colnames(Y) <- c("trait1", "trait2")

approx <- lm_phen_cov_laplace(Y, X, draws = 2000, seed = 41)
sampled <- lm_phen_cov(Y, X, iter = 1000, warmup = 500, chains = 2,
                       cores = 2, seed = 41, refresh = 0)

print(approx)
print(sampled$summary(c("P", "corrP"))[, c("variable", "mean", "sd", "rhat")])
print(sampled$diagnostic_summary())
cat("\nThe approximation is Gaussian on the unconstrained scale, not a",
    "posterior sample. Compare estimates and check the sampler diagnostics;",
    "agreement in this example does not guarantee agreement in other models.",
    "Use the sampler for uncertainty and inference.\n")
