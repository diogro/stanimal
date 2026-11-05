# First look at the Calomys skull data: unstructured against factor
#
# The same eight traits as the mice, four face and four neurocranium, so that a
# difference between the two datasets is a difference of design and not of
# traits. Eight traits also keep the unstructured G feasible at 36 covariance
# parameters, so the two models are compared on equal terms.
#
# The slab is 1, which is now the package default: through the row normalisation
# it states that the prior expects a communality of 0.5, against 0.8 at a slab of
# 2. With mean |r| = 0.667 here the data will probably want more than the prior
# allows, and that is worth seeing rather than hiding behind a permissive slab.
#
# Run from the repository root:
#   Rscript data-raw/calomys-fit.R

suppressMessages(pkgload::load_all("package", quiet = TRUE))

d <- readRDS("data-raw/calomys.rds")
skull <- d$skull
meta <- c("id", "sire", "dam", "sex", "age", "litter", "experiment")
traits <- c("NSL_ZI", "PM_MT", "NSL_ZS", "IS_PNS",
            "PT_APET", "BR_APET", "PNS_APET", "OPI_LD")
stopifnot(all(traits %in% names(skull)))

skull <- skull[, c(meta, traits)]
skull$log_age <- log(skull$age)
skull$experiment <- factor(skull$experiment)
rownames(skull) <- skull$id

design <- animal_design(
  stats::reformulate(c("sex", "log_age", "experiment"),
                     response = sprintf("cbind(%s)", paste(traits, collapse = ", "))),
  data = skull, pedigree = d$skull_pedigree)

cat("Y:", nrow(design$Y), "x", ncol(design$Y), " X:", ncol(design$X),
    " A:", nrow(design$A), "\n")

ess_min <- function(fit) {
  min(posterior::summarize_draws(fit$draws(variables = "G"),
                                 ess_bulk = posterior::ess_bulk)$ess_bulk,
      na.rm = TRUE)
}

fit_one <- function(label, ...) {
  elapsed <- system.time(
    fit <- lmm_animal(design$Y, design$X, design$A, iter = 1000, warmup = 500,
                      chains = 4, cores = 4, refresh = 0, seed = 20261105L, ...)
  )[["elapsed"]]
  diag_s <- fit$diagnostic_summary(quiet = TRUE)
  h2 <- draws_mean(fit, "h2", traits = traits)
  corrG <- draws_mean(fit, "corrG", traits = traits)
  comm <- tryCatch(draws_mean(fit, "communality", traits = traits),
                   error = function(e) NULL)
  cat(sprintf("\n== %s ==\n", label))
  cat(sprintf("h2 %.3f [%.3f, %.3f]  mean|corrG| %.3f  %s\n",
              mean(h2), min(h2), max(h2),
              mean(abs(corrG[lower.tri(corrG)])),
              if (is.null(comm)) "" else
                sprintf("communality %.3f", mean(comm))))
  cat(sprintf("ESS %.0f  Rhat %.3f  treedepth %d  divergences %d  %.0f s\n",
              ess_min(fit), max(fit$summary()$rhat, na.rm = TRUE),
              sum(diag_s$num_max_treedepth), sum(diag_s$num_divergent), elapsed))
  saveRDS(fit, sprintf("data-raw/calomys8-%s.rds", label))
  invisible(fit)
}

fit_one("unstructured")
fit_one("factor-M2-slab1", genetic_structure = "factor", n_factors = 2,
        hs_slab_scale = 1)
