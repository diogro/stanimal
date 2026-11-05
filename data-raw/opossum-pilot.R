# Opossum pilot: eight traits chosen for a spread of phenotypic correlations
#
# The point of this pilot is the trait selection. Both earlier real datasets
# failed for the same reason: every pair of trait correlations sat near 0.9, so
# the covariance had one dominant direction, a second factor had nothing to
# explain, and the decomposition saturated. Here the eight are chosen so that the
# 28 pairwise correlations span a wide range and the first eigenvalue is not
# dominant, which gives a factor model something to find and something to fail at.
#
# The relationship matrix is genomic, built from the 3,696 coded markers, because
# the shipped pedigree kinship has every pair related at about 0.77 and cannot
# separate additive variance from a shared colony effect.
#
# The traits are used as measured: Porto et al. analysed these distances
# untransformed, and the published P matrix confirms it (correlation 0.834 with
# the raw variances against 0.046 with the logged ones).
#
# Run from the repository root:
#   Rscript data-raw/opossum-pilot.R

suppressMessages(pkgload::load_all("package", quiet = TRUE))

src <- "/Users/diogro/Downloads/4055961"

## Data -----------------------------------------------------------------------

cranial <- read.csv(file.path(src, "phenotype_covariates", "cranial_data.csv"),
                    stringsAsFactors = FALSE)

landmarks <- c("IS", "NSL", "PNS", "PM", "ZS", "ZI", "MT", "NA", "BR", "PT",
               "APET", "BA", "EAM", "ZYGO", "TSP", "LD", "AS", "OPI", "JP")
landmarks <- landmarks[order(-nchar(landmarks))]
porto_name <- function(x) {
  vapply(x, function(one) {
    one <- sub("[0-9]+$", "", one)
    for (left in landmarks)
      if (startsWith(one, left) && nchar(one) > nchar(left))
        return(paste0(left, "_", substring(one, nchar(left) + 1L)))
    one
  }, character(1), USE.NAMES = FALSE)
}
measured <- names(cranial)[4:38]
traits_all <- porto_name(measured)

# The marker files carry no animal IDs, so animals are matched by row order. The
# cranial file has 553 rows and each marker file has 553, which is the only
# check available.
unzip(file.path(src, "coded_marker_data.zip"), exdir = tempdir())
marker_files <- list.files(file.path(tempdir(), "coded_marker_data"),
                           pattern = "SAS_ADD", full.names = TRUE)
stopifnot(length(marker_files) == 8L)
markers <- do.call(cbind, lapply(marker_files, function(f)
  as.matrix(read.csv(f, header = FALSE))))
stopifnot(nrow(markers) == nrow(cranial))

freq <- (colMeans(markers) + 1) / 2
Z <- sweep(markers, 2, 2 * freq - 1, "-")
GRM <- tcrossprod(Z) / ncol(markers)
dimnames(GRM) <- list(cranial$Animal, cranial$Animal)

complete <- stats::complete.cases(cranial[measured]) &
  !is.na(cranial$AgeBI) & !is.na(cranial$SexBI)
cat("animals:", nrow(cranial), " complete records:", sum(complete),
    " traits:", length(traits_all), "\n")

## Residual phenotypic correlations after the covariates ----------------------

dat <- cranial[complete, c("Animal", "AgeBI", "SexBI")]
dat[traits_all] <- as.matrix(cranial[complete, measured])
resid <- sapply(traits_all, function(trait)
  stats::residuals(stats::lm(stats::reformulate(c("AgeBI", "SexBI"), trait),
                             data = dat)))
R <- stats::cor(resid)
dimnames(R) <- list(traits_all, traits_all)

## Choose eight for a wide spread of correlations -----------------------------

score_set <- function(set) {
  sub <- R[set, set][lower.tri(diag(length(set)))]
  ev <- eigen(R[set, set], symmetric = TRUE, only.values = TRUE)$values
  c(spread = stats::sd(abs(sub)),
    strong = mean(abs(sub) > 0.5),
    weak = mean(abs(sub) < 0.2),
    pc1 = ev[1] / sum(ev))
}

set.seed(20261105)
best <- NULL; best_score <- -Inf
for (i in 1:50000) {
  set <- sample(traits_all, 8L)
  s <- score_set(set)
  # A set earns its place by spreading the correlations, keeping a share of weak
  # ones, and staying away from a single dominant axis.
  value <- s["spread"] + 0.5 * s["weak"] - 2 * max(0, s["pc1"] - 0.60)
  if (value > best_score) { best_score <- value; best <- set; best_s <- s }
}
traits <- sort(best)
cat("\nSelected eight:", paste(traits, collapse = ", "), "\n")
cat(sprintf("spread sd|r| %.3f  share |r|>0.5 %.2f  share |r|<0.2 %.2f  PC1 %.2f\n",
            best_s["spread"], best_s["strong"], best_s["weak"], best_s["pc1"]))
print(round(R[traits, traits], 2))
cat("off-diagonal |r| quantiles:",
    paste(round(quantile(abs(R[traits, traits][lower.tri(diag(8))]),
                         c(0, .25, .5, .75, 1)), 2), collapse = " "), "\n")

## Fit -----------------------------------------------------------------------

pilot <- data.frame(id = dat$Animal, sex = factor(dat$SexBI),
                    age = dat$AgeBI, stringsAsFactors = FALSE)
pilot[traits] <- dat[traits]
rownames(pilot) <- pilot$id
A <- GRM[pilot$id, pilot$id]

# The design matrix and response are built directly, since the relationship
# matrix is genomic and animal_design() would want a pedigree to build one.
Y <- as.matrix(pilot[traits])
rownames(Y) <- pilot$id
X <- stats::model.matrix(~ sex + age, data = pilot)
A <- GRM[pilot$id, pilot$id]

ess_min <- function(fit) min(posterior::summarize_draws(
  fit$draws(variables = "G"), ess_bulk = posterior::ess_bulk)$ess_bulk, na.rm = TRUE)

fit_one <- function(label, ...) {
  elapsed <- system.time(
    fit <- lmm_animal(Y, X, A, iter = 1000, warmup = 500, chains = 4, cores = 4,
                      refresh = 0, seed = 20261105L, ...)
  )[["elapsed"]]
  d <- fit$diagnostic_summary(quiet = TRUE)
  h2 <- draws_mean(fit, "h2", traits = traits)
  corrG <- draws_mean(fit, "corrG", traits = traits)
  comm <- tryCatch(draws_mean(fit, "communality", traits = traits),
                   error = function(e) NULL)
  cat(sprintf("\n== %s ==\nh2 %.3f [%.3f, %.3f]  mean|corrG| %.3f  %s\n",
              label, mean(h2), min(h2), max(h2),
              mean(abs(corrG[lower.tri(corrG)])),
              if (is.null(comm)) "" else sprintf("communality %.3f", mean(comm))))
  cat(sprintf("ESS %.0f  Rhat %.3f  treedepth %d  divergences %d  %.0f s\n",
              ess_min(fit), max(fit$summary()$rhat, na.rm = TRUE),
              sum(d$num_max_treedepth), sum(d$num_divergent), elapsed))
  saveRDS(fit, sprintf("data-raw/opossum-pilot-%s.rds", label))
  invisible(fit)
}

fit_one("unstructured")
fit_one("factor-M2", genetic_structure = "factor", n_factors = 2, hs_slab_scale = 1)
