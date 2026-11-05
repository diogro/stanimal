# Prototype: the factor model with the specific genetic variance free
#
# Compares the shipped factor model (loadings normalised to a correlation
# matrix, so the specific genetic variance is 1 / (1 + ||raw||^2) and the trait
# heritability is a parameter) with animalModelGFactorPsi.stan, where
# G = Lambda Lambda' + diag(psi) and the heritability is derived. See Runcie and
# Mukherjee (2013) for the parameterisation of G, with unstructured rather than
# diagonal residuals here.
#
# The point of the comparison is stability: does the reported heritability still
# move when the slab scale moves, and does it still move when the trait selection
# changes from six facial and two neurocranial distances to a balanced four and
# four (Porto et al. 2009, Table 2)?
#
# Run from the repository root:
#   Rscript data-raw/prototype-factor-psi.R
# Writes data-raw/prototype-factor-psi.csv.

suppressMessages(pkgload::load_all("package", quiet = TRUE))

data(ratones_morpho)

# ratones_morpho records each mouse's dam and sire, so the pedigree can be
# rebuilt without a separate file; parents that were never measured are added
# as founders by makeA().
ratones_pedigree <- data.frame(id = as.character(ratones_morpho$ID),
                               sire = as.character(ratones_morpho$PADRE),
                               dam = as.character(ratones_morpho$MADRE),
                               stringsAsFactors = FALSE)

selections <- list(
  "6 face : 2 neuro" = c("NSL_ZI", "PM_MT", "IS_PNS", "NSL_ZS", "PT_EAM",
                         "PT_ZYGO", "PT_APET", "PM_ZI"),
  "4 + 4" = c("NSL_ZI", "PM_MT", "NSL_ZS", "IS_PNS",
              "PT_APET", "BR_APET", "PNS_APET", "OPI_LD"))

psi_model <- normalizePath("package/inst/stan/animalModelGFactorPsi.stan")

build <- function(traits) {
  dat <- ratones_morpho
  dat[traits] <- log(dat[traits])
  design <- suppressWarnings(animal_design(
    stats::reformulate(c("SEX", "AGE", "LIN"),
                       response = sprintf("cbind(%s)",
                                          paste(traits, collapse = ", "))),
    data = dat, pedigree = ratones_pedigree))
  resid <- sapply(traits, function(trait) {
    stats::residuals(stats::lm(
      stats::reformulate(c("SEX", "AGE", "LIN"), trait),
      data = cbind(dat[traits], dat[c("SEX", "AGE", "LIN")])))
  })
  Rp <- stats::cor(resid)
  dimnames(Rp) <- list(traits, traits)
  list(Y = design$Y, X = design$X, A = design$A, traits = traits, corrP = Rp)
}

ess_min <- function(fit) {
  min(posterior::summarize_draws(fit$draws(variables = "G"),
                                 ess_bulk = posterior::ess_bulk)$ess_bulk,
      na.rm = TRUE)
}

fit_one <- function(label, d, ...) {
  elapsed <- system.time(
    fit <- lmm_animal(d$Y, d$X, d$A, iter = 2000, warmup = 1000,
                      chains = 4, cores = 4, refresh = 0, seed = 20261105L, ...)
  )[["elapsed"]]
  diag_s <- fit$diagnostic_summary(quiet = TRUE)
  traits <- d$traits
  G <- draws_mean(fit, "G", traits = traits)
  E <- draws_mean(fit, "E", traits = traits)
  h2 <- draws_mean(fit, "h2", traits = traits)
  corrG <- draws_mean(fit, "corrG", traits = traits)
  comm <- tryCatch(draws_mean(fit, "communality", traits = traits),
                   error = function(e) NULL)

  # The heritability is the weight on the genetic side of
  # corrP_ij = sqrt(h2_i h2_j) corrG_ij + sqrt((1-h2_i)(1-h2_j)) corrE_ij,
  # so the share of the observed phenotypic correlation attributed to G is the
  # quantity that has to be stable if the heritability is to mean anything.
  sE <- sqrt(diagonal_of(E))
  sG <- sqrt(diagonal_of(G))
  corrE <- E / outer(sE, sE)
  wG <- outer(h2, h2, function(a, b) sqrt(a * b))
  wE <- outer(1 - h2, 1 - h2, function(a, b) sqrt(a * b))
  impl <- wG * corrG + wE * corrE
  lower <- lower.tri(d$corrP)

  row <- data.frame(
    fit = label,
    h2_mean = mean(h2), h2_min = min(h2), h2_max = max(h2),
    h2_sd_between_traits = stats::sd(h2),
    communality_mean = if (is.null(comm)) NA_real_ else mean(comm),
    corrG_mean_abs = mean(abs(corrG[lower])),
    share_corr_in_G = sum(wG[lower] * corrG[lower]) / sum(d$corrP[lower]),
    corrP_error = max(abs(impl[lower] - d$corrP[lower])),
    min_ess = ess_min(fit),
    max_rhat = max(fit$summary()$rhat, na.rm = TRUE),
    treedepth = sum(diag_s$num_max_treedepth),
    divergences = sum(diag_s$num_divergent),
    time_s = elapsed)
  cat("\n"); print(row, row.names = FALSE, digits = 3)
  row
}

diagonal_of <- function(M) diag(M)

rows <- list()
for (design in names(selections)) {
  d <- build(selections[[design]])
  cat("\n====", design, "====\n")
  rows[[length(rows) + 1L]] <- fit_one(
    paste(design, "| current, slab 1"), d,
    genetic_structure = "factor", n_factors = 2, hs_slab_scale = 1)
  rows[[length(rows) + 1L]] <- fit_one(
    paste(design, "| current, slab 2"), d,
    genetic_structure = "factor", n_factors = 2, hs_slab_scale = 2)
  rows[[length(rows) + 1L]] <- fit_one(
    paste(design, "| psi, slab 1"), d,
    genetic_structure = "factor", n_factors = 2, hs_slab_scale = 1,
    model_file = psi_model)
  rows[[length(rows) + 1L]] <- fit_one(
    paste(design, "| psi, slab 2"), d,
    genetic_structure = "factor", n_factors = 2, hs_slab_scale = 2,
    model_file = psi_model)
  rows[[length(rows) + 1L]] <- fit_one(paste(design, "| unstructured"), d)
  write.csv(do.call(rbind, rows), "data-raw/prototype-factor-psi.csv",
            row.names = FALSE)
}
cat("\nFull table\n")
print(do.call(rbind, rows), row.names = FALSE, digits = 3)
