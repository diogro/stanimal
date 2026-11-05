# Opossum pilot, second pass: the h2 parametrisation and one to four factors
#
# The first pass gave h2 of 0.195 (unstructured) and 0.241 (two factors), with
# overlapping intervals and clean convergence, on eight traits chosen for a wide
# spread of phenotypic correlations. This pass asks two follow-ups: re-read the
# unstructured heritability under the h2 parametrisation, since one trait of
# eight sat near the boundary under the variance one, and fit one to four factors
# to see whether the correlated cluster really is a single direction.
#
# Run from the repository root:
#   Rscript data-raw/opossum-models.R

suppressMessages(pkgload::load_all("package", quiet = TRUE))

src <- "/Users/diogro/Downloads/4055961"
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
cranial[porto_name(measured)] <- cranial[measured]

# Same eight as the pilot, chosen there for the spread of correlations.
traits <- c("BA_OPI", "BR_PT", "PT_APET", "PT_BA", "PT_EAM", "ZI_MT", "ZI_TSP",
            "ZS_ZI")

unzip(file.path(src, "coded_marker_data.zip"), exdir = tempdir())
marker_files <- list.files(file.path(tempdir(), "coded_marker_data"),
                           pattern = "SAS_ADD", full.names = TRUE)
markers <- do.call(cbind, lapply(marker_files, function(f)
  as.matrix(read.csv(f, header = FALSE))))
freq <- (colMeans(markers) + 1) / 2
GRM <- tcrossprod(sweep(markers, 2, 2 * freq - 1, "-")) / ncol(markers)
dimnames(GRM) <- list(cranial$Animal, cranial$Animal)

ok <- stats::complete.cases(cranial[traits]) &
  !is.na(cranial$AgeBI) & !is.na(cranial$SexBI)
d <- data.frame(id = cranial$Animal[ok], sex = factor(cranial$SexBI[ok]),
                age = cranial$AgeBI[ok], stringsAsFactors = FALSE)
d[traits] <- cranial[ok, traits]
Y <- as.matrix(d[traits]); rownames(Y) <- d$id
X <- stats::model.matrix(~ sex + age, data = d)
A <- GRM[d$id, d$id]
cat("Y:", nrow(Y), "x", ncol(Y), " complete records\n")

ess_min <- function(fit) min(posterior::summarize_draws(
  fit$draws(variables = "G"), ess_bulk = posterior::ess_bulk)$ess_bulk, na.rm = TRUE)

fit_one <- function(label, ...) {
  elapsed <- system.time(
    fit <- lmm_animal(Y, X, A, iter = 1000, warmup = 500, chains = 4, cores = 4,
                      refresh = 0, seed = 20261105L, ...)
  )[["elapsed"]]
  dg <- fit$diagnostic_summary(quiet = TRUE)
  h2 <- draws_mean(fit, "h2", traits = traits)
  corrG <- draws_mean(fit, "corrG", traits = traits)
  comm <- tryCatch(draws_mean(fit, "communality", traits = traits),
                   error = function(e) NULL)
  L <- tryCatch(draws_mean(fit, "loadings", traits = traits), error = function(e) NULL)
  elpd <- tryCatch({
    l <- suppressWarnings(fit_loo(fit))
    c(mean(l$estimates["elpd_loo", "Estimate"]),
      mean(l$diagnostics$pareto_k > 0.7))
  }, error = function(e) c(NA, NA))

  cat(sprintf("\n== %s ==\nh2 %.3f [%.3f, %.3f]  mean|corrG| %.3f  %s elpd %.1f  k>0.7 %d\n",
              label, mean(h2), min(h2), max(h2),
              mean(abs(corrG[lower.tri(corrG)])),
              if (is.null(comm)) "" else sprintf("comm %.3f ", mean(comm)),
              elpd[1], as.integer(elpd[2])))
  cat(sprintf("ESS %.0f  Rhat %.3f  treedepth %d  divergences %d  %.0f s\n",
              ess_min(fit), max(fit$summary()$rhat, na.rm = TRUE),
              sum(dg$num_max_treedepth), sum(dg$num_divergent), elapsed))
  if (!is.null(L)) {
    cat("loadings:\n")
    print(round(as.matrix(L), 2))
  }
  saveRDS(fit, sprintf("data-raw/opossum-%s.rds", label))
  invisible(list(fit = fit, elpd = elpd))
}

res <- list()
res[["us-h2"]] <- fit_one("us-h2", parametrization = "h2")
for (M in 1:4)
  res[[paste0("factor-M", M)]] <- fit_one(paste0("factor-M", M),
                                         genetic_structure = "factor",
                                         n_factors = M, hs_slab_scale = 1)

cat("\nelpd summary\n")
print(round(sapply(res, function(x) x$elpd[1]), 1))
