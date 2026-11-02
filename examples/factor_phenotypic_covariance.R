# Sparse factor phenotypic covariance on the full mouse morphometric data
#
# The 35 log-transformed mouse traits are dominated by a single size gradient:
# PC1 alone takes about 40% of the variance and the eigenvalues then fall away
# quickly. That is the low-rank structure a factor model describes, so
# lm_phen_cov() estimates P as a few sparse factors plus trait-specific
# variances instead of a full 35 x 35 correlation matrix.
#
# The factors stay on the observed log-trait scale on purpose. Rotating the
# traits first (size plus shape scores, see size_shape()) spreads every loading
# over all 35 traits, which destroys the sparsity the prior is there to find,
# and it moves the low-rank structure out of P and into the rotation itself.
#
# Run from the project root:
#   Rscript examples/factor_phenotypic_covariance.R

if (!requireNamespace("cmdstanr", quietly = TRUE))
  stop("Install cmdstanr and CmdStan first: https://mc-stan.org/cmdstanr/")

library(stanimal)

ratonesdf <- ratones_morpho
traits <- names(ratonesdf)[12:46]
Y <- log(as.matrix(ratonesdf[, traits]))
colnames(Y) <- traits

cat(nrow(Y), "individuals,", ncol(Y), "log traits\n")
cat("PC1 variance share:", round(pc1_variance_share(Y), 3), "\n")
eigenvalues <- eigen(stats::cor(Y), symmetric = TRUE, only.values = TRUE)$values
cat("first five eigenvalues:", round(eigenvalues[1:5], 2), "\n")

# The sex, age and selection-line effects are the fixed effects; there is no
# relatedness matrix in this example.
X <- stats::model.matrix(~ SEX + AGE + LIN, data = ratonesdf)

# Three factors: the size axis plus two shape contrasts, which is where the
# eigenvalues stop standing out. Extra factors only re-split variance the
# existing ones already carry, and the rotation they introduce makes the factor
# summaries multimodal without improving P. The mild column decay lets the broad
# first factor come before the narrower ones without forcing an order.
fit <- lm_phen_cov(
  Y, X,
  covariance_structure = "factor",
  n_factors = 3,
  hs_global_scale = 2, hs_slab_scale = 2,
  column_decay = 0.8,
  iter = 2000, warmup = 1000, chains = 4, cores = 4,
  refresh = 0, show_messages = FALSE
)

# P is already on the log-trait scale, so it is directly comparable with the
# OLS residual covariance, the unpenalised estimate of the same quantity. The
# error is reported relative to the size of the covariance: these variances are
# around 0.001, so an absolute RMSE of a few thousandths is a large relative
# error and says nothing on its own.
P_hat <- draws_mean(fit, "P")
P_ols <- stats::cov(stats::lm.fit(X, Y)$residuals)
cat("P relative error against the OLS residual covariance:",
    round(sqrt(mean((P_hat - P_ols)^2)) / sqrt(mean(P_ols^2)), 3), "\n")
cat("first five variances:\n")
print(round(rbind(factor = diag(P_hat)[1:5], ols = diag(P_ols)[1:5]), 4))

# The horseshoe shrinks, so the fitted correlations should sit slightly closer
# to zero than the OLS ones. A large gap the other way would mean the prior is
# inflating the common structure instead of smoothing it.
corr_hat <- stats::cov2cor(P_hat)
corr_ols <- stats::cov2cor(P_ols)
cat("mean |r|:", round(mean(abs(corr_hat[upper.tri(corr_hat)])), 3),
    "for the factor model and",
    round(mean(abs(corr_ols[upper.tri(corr_ols)])), 3), "for OLS\n")

# The factor summaries only mean something because the factors live on the
# observed trait scale: communality is the share of a trait's variance taken by
# the common factors, uniqueness is its complement, and a small uniqueness for a
# handful of traits is the sparse module the prior is looking for. Loadings and
# factor magnitudes are not rotation invariant, so read them only after P has
# converged; communality and uniqueness are invariant.
lambda <- draws_mean(fit, "loadings")
cat("\nposterior mean loadings, first five traits:\n")
print(round(lambda[1:5, ], 3))
# factor_magnitude is not rotation invariant. These three factors are of
# comparable size, so their labels exchange between chains and the per-factor
# magnitudes mix the modes: that is what the high Rhat and small ESS below
# report, and it is not a sign that P failed. The rotation-invariant summary of
# the same information is the total common variance, sum(communality).
cat("\ntotal common variance (sum of communality):",
    round(sum(draws_mean(fit, "communality")), 2), "\n")
print(fit$summary("factor_magnitude")[, c("variable", "mean", "sd",
                                           "rhat", "ess_bulk")])
print(fit$summary("communality")[1:5, c("mean", "sd", "rhat", "ess_bulk")])
print(fit$summary("uniqueness")[1:5, c("mean", "sd", "rhat", "ess_bulk")])
print(fit$diagnostic_summary())

# Heatmaps: the fitted correlations, the OLS ones and the loadings that
# produced them.
dir.create("examples/figures", recursive = TRUE, showWarnings = FALSE)
pdf("examples/figures/factor_phenotypic_covariance.pdf",
    width = 12, height = 12)
print(plot_corr(corr_hat, "Factor P: posterior mean correlations"))
print(plot_corr(corr_ols, "OLS residual correlations"))
print(plot_corr(t(lambda), "Posterior mean loadings (correlation scale)"))
dev.off()
