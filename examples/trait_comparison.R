# The same animal model fitted with several engines, for a chosen number of
# traits
#
# Simulates K genetically correlated traits with a known G and E on a
# five-generation pedigree, fits the same multivariate animal model with
# stanimal, MCMCglmm, BGLR, sommer, lme4breeding and ASReml-R, and compares the
# estimates with the simulated truth and the time each fit took.
#
# K is the first command-line argument (default 5):
#   Rscript examples/trait_comparison.R       # five traits
#   Rscript examples/trait_comparison.R 3     # three traits
#
# stanimal, MCMCglmm and BGLR sample the posterior; sommer, lme4breeding and
# ASReml-R maximize the REML likelihood. The open-source engines are installed
# below if they are missing; ASReml-R is licensed and has to be installed
# separately.

# Install the engines that are missing (stanimal is assumed installed)
pkgs <- c("MCMCglmm", "BGLR", "sommer", "lme4breeding")
new <- setdiff(pkgs, rownames(utils::installed.packages()))
if (length(new))
  install.packages(new, repos = "https://cloud.r-project.org")

# ASReml-R is licensed and not on CRAN: only check whether it is available,
# never try to install it
has_asreml <- requireNamespace("asreml", quietly = TRUE)

library(stanimal)

## Number of traits ----------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
n_traits <- if (length(args)) as.integer(args[1L]) else 5L
if (is.na(n_traits) || n_traits < 2L)
  stop("the number of traits must be a whole number of at least 2.",
       call. = FALSE)
traits <- paste0("t", seq_len(n_traits))

## Simulation ----------------------------------------------------------------

set.seed(1234)

n_gen <- 5; n_breed <- 100
ped <- data.frame(id = seq_len(n_gen * n_breed),
                  sire = NA_integer_, dam = NA_integer_,
                  sex = NA_character_, stringsAsFactors = FALSE)
# Sex-consistent parents: a male can only be a sire and a female only a dam.
# This also keeps the pedigree less regular, so REML is better conditioned.
ped$sex[seq_len(n_breed)] <- sample(rep(c("M", "F"), each = n_breed / 2))
for (g in 2:n_gen) {
  kids <- ((g - 1) * n_breed + 1):(g * n_breed)
  pool <- ped[ped$id < min(kids), ]
  ped$sire[kids] <- sample(pool$id[pool$sex == "M"], n_breed, TRUE)
  ped$dam[kids]  <- sample(pool$id[pool$sex == "F"], n_breed, TRUE)
  ped$sex[kids]  <- sample(c("M", "F"), n_breed, TRUE)
}
A <- makeA(ped)

# Simulated truth from random_covariance(): the default LKJ construction shrinks
# correlations toward zero, so G and E are well conditioned and REML does not
# sit on the boundary.
set.seed(4321)
G <- random_covariance(n_traits, min_var = 0.7, max_var = 1.3,
                       trait_names = traits)
E <- random_covariance(n_traits, min_var = 0.7, max_var = 1.3,
                       trait_names = traits)

n <- nrow(A)
X <- cbind(Intercept = 1, sex = rbinom(n, 1, 0.5), age = round(rnorm(n, 3, 0.5), 2))
B <- rbind(rep(10, n_traits), rep(1.5, n_traits), rep(0.4, n_traits))

Y <- X %*% B +
  t(chol(A)) %*% matrix(rnorm(n * n_traits), n) %*% chol(G) +
  matrix(rnorm(n * n_traits), n) %*% chol(E)
colnames(Y) <- traits

dat <- data.frame(id = factor(ped$id), sex = X[, "sex"], age = X[, "age"], Y)
# animal_design() needs the individual IDs as row names to line the pedigree up
rownames(dat) <- as.character(dat$id)

## Comparison helpers --------------------------------------------------------

# ASReml reports us() parameters as the lower triangle, row by row
lower2mat <- function(x, k) {
  m <- matrix(0, k, k)
  rows <- rep(seq_len(k), seq_len(k))
  cols <- sequence(seq_len(k))
  m[cbind(rows, cols)] <- x
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  m
}

corr_of <- function(S) {
  if (nrow(S) < 2L) return(numeric(0))
  cov2cor(S)[upper.tri(S)]
}

# ASReml's vm() wants a sparse inverse as a three-column matrix: row, column
# and value for the lower triangle, sorted row-major. Building it from A avoids
# asreml::ainverse(), which fails intermittently on some pedigrees.
sparse_ginv <- function(A) {
  inverse <- solve(A)
  inverse[abs(inverse) < 1e-10 * max(abs(inverse))] <- 0
  lower <- lower.tri(inverse, diag = TRUE) & inverse != 0
  index <- which(lower, arr.ind = TRUE)
  index <- index[order(index[, 1], index[, 2]), , drop = FALSE]
  out <- cbind(Row = index[, 1], Column = index[, 2],
               Ainverse = inverse[index])
  attr(out, "rowNames") <- rownames(A)
  attr(out, "geneticGroups") <- c(0, 0)
  class(out) <- c(class(out), "ginv")
  out
}

h2_true <- diag(G) / diag(G + E)
rG_true <- corr_of(G)

# Relative MSE divides the elementwise mean squared error by the mean squared
# truth on the same scale. That is what makes the columns comparable across
# quantities: a covariance entry can be 0.001, where a small absolute error is a
# large relative one, while a heritability near 1 has the opposite problem.
relative_mse <- function(estimate, truth)
  mean((estimate - truth)^2) / mean(truth^2)

# One row per engine: relative MSE of G, E, the heritabilities and the genetic
# correlations, the minimum ESS of G, ESS per retained draw, and time.
compare <- function(engine, G_hat, E_hat, elapsed = NA_real_,
                    min_ess = NA_real_, n_draws = NA_real_) {
  ok <- !is.null(G_hat) && !is.null(E_hat) &&
    !anyNA(G_hat) && !anyNA(E_hat)
  if (!ok)
    return(data.frame(engine = engine, G_relative_mse = NA_real_,
                      E_relative_mse = NA_real_,
                      h2_relative_mse = NA_real_,
                      corrG_relative_mse = NA_real_,
                      min_ess = min_ess, ess_ratio = min_ess / n_draws,
                      time_s = elapsed))
  data.frame(engine = engine,
             G_relative_mse = relative_mse(G_hat, G),
             E_relative_mse = relative_mse(E_hat, E),
             h2_relative_mse = relative_mse(
               diag(G_hat) / diag(G_hat + E_hat), h2_true),
             corrG_relative_mse = relative_mse(corr_of(G_hat), rG_true),
             min_ess = min_ess, ess_ratio = min_ess / n_draws,
             time_s = elapsed)
}

# Smallest bulk ESS across the unique elements of G. The draws matrix has one
# row per retained draw and one column per parameter.
ess_min <- function(draws) {
  draws <- as.matrix(draws)
  if (!nrow(draws) || !ncol(draws)) return(NA_real_)
  values <- vapply(seq_len(ncol(draws)),
                   function(j) posterior::ess_bulk(draws[, j]), numeric(1))
  if (all(is.na(values))) return(NA_real_)
  min(values, na.rm = TRUE)
}

results <- list(compare("truth", G, E))

# One row per engine, filled in as each fit runs, so an engine that is not
# installed (or whose fit is skipped) still appears in the table, with NA
# estimates and no time.
engines <- c("stanimal", "MCMCglmm", "BGLR", "sommer", "lme4breeding",
             "asreml")
na_mat <- matrix(NA_real_, n_traits, n_traits)
for (engine in engines)
  results[[engine]] <- compare(engine, na_mat, na_mat)

## stanimal: Bayesian, breeding values marginalized out of the likelihood ----

response <- sprintf("cbind(%s)", paste(traits, collapse = ", "))
design <- animal_design(reformulate(c("sex", "age"), response = response),
                        dat, pedigree = ped)
time_stan <- system.time(
  fit_stan <- lmm_animal(design$Y, design$X, design$A,
                         iter = 2000, warmup = 1000, chains = 4, cores = 4,
                         refresh = 0))
g_pairs <- which(lower.tri(diag(n_traits), diag = TRUE), arr.ind = TRUE)
g_draws <- fit_draws(fit_stan)$G
g_matrix <- vapply(seq_len(nrow(g_pairs)), function(i)
  g_draws[g_pairs[i, 1], g_pairs[i, 2], ], numeric(dim(g_draws)[3]))
results$stanimal <- compare("stanimal", draws_mean(fit_stan, "G"),
                            draws_mean(fit_stan, "E"),
                            elapsed = time_stan[["elapsed"]],
                            min_ess = ess_min(g_matrix),
                            n_draws = nrow(g_matrix))

## MCMCglmm: Bayesian, pedigree supplied as the inverse of A ----------------

library(MCMCglmm)
prior <- list(G = list(G = list(V = diag(n_traits) * 0.5, nu = n_traits)),
              R = list(R = list(V = diag(n_traits) * 0.5, nu = n_traits)))
time_mcmc <- system.time(
  fit_mcmc <- MCMCglmm(as.formula(paste(response,
                                        "~ trait + trait:sex + trait:age - 1")),
                       random   = ~ us(trait):id,
                       rcov     = ~ us(trait):units,
                       ginverse = list(id = MCMCglmm::inverseA(ped[, 1:3])$Ainv),
                       family   = rep("gaussian", n_traits),
                       prior = prior, data = dat,
                       nitt = 13000, burnin = 3000, thin = 10, verbose = FALSE))
# VCV holds the K x K G matrix in column order, then the K x K residual matrix
k2 <- n_traits^2
mcmc_draws <- as.matrix(fit_mcmc$VCV)[, seq_len(k2), drop = FALSE]
results$MCMCglmm <- compare(
  "MCMCglmm",
  matrix(colMeans(fit_mcmc$VCV[, seq_len(k2)]), n_traits),
  matrix(colMeans(fit_mcmc$VCV[, k2 + seq_len(k2)]), n_traits),
  elapsed = time_mcmc[["elapsed"]],
  min_ess = ess_min(mcmc_draws), n_draws = nrow(mcmc_draws))

## BGLR: Bayesian multitrait, pedigree relationship as an RKHS kernel -------

if (requireNamespace("BGLR", quietly = TRUE)) {
  X_bglr <- model.matrix(~ sex + age, dat)
  A_dat <- A[match(as.character(dat$id), rownames(A)),
             match(as.character(dat$id), rownames(A))]
  # Weak inverse-Wishart priors with prior mean I, for G and for the residual
  prior_cov <- list(type = "UN", df0 = n_traits + 3,
                    S0 = 2 * diag(n_traits))
  bglr_prefix <- file.path(tempdir(), "bglr_traits_")
  time_bglr <- system.time(
    fit_bglr <- BGLR::Multitrait(
      y = as.matrix(dat[traits]),
      ETA = list(list(X = X_bglr, model = "FIXED"),
                 list(K = A_dat, model = "RKHS", Cov = prior_cov)),
      intercept = FALSE, resCov = prior_cov,
      nIter = 13000, burnIn = 3000, thin = 10,
      saveAt = bglr_prefix, verbose = FALSE))
  bglr_omega <- Sys.glob(paste0(bglr_prefix, "Omega_*.dat"))
  bglr_ess <- NA_real_
  bglr_draws <- NA_real_
  if (length(bglr_omega)) {
    bglr_samples <- as.matrix(read.table(bglr_omega[1]))
    bglr_ess <- ess_min(bglr_samples)
    bglr_draws <- nrow(bglr_samples)
  }
  results$BGLR <- compare("BGLR", fit_bglr$ETA[[2]]$Cov$Omega,
                          fit_bglr$resCov$R,
                          elapsed = time_bglr[["elapsed"]],
                          min_ess = bglr_ess, n_draws = bglr_draws)
} else {
  message("BGLR is not installed; install.packages(\"BGLR\") to add it")
}

## sommer: REML through mmes(), in long format -------------------------------

library(sommer)
# mmes() takes the inverse relationship matrix, marked as such, and a
# long-format response: usm(trait) estimates the trait covariance and ism(id)
# applies A^-1 to the animal effect.
Ainv <- solve(A)
attr(Ainv, "inverse") <- TRUE
long_sommer <- data.frame(
  id = factor(rep(dat$id, n_traits)),
  trait = factor(rep(traits, each = nrow(dat)), levels = traits),
  value = unlist(dat[traits], use.names = FALSE),
  sex = rep(dat$sex, n_traits),
  age = rep(dat$age, n_traits),
  record = rep(seq_len(nrow(dat)), n_traits))
time_sommer <- system.time(
  fit_sommer <- mmes(
    value ~ 0 + trait + trait:sex + trait:age,
    random = ~ vsm(usm(trait), ism(id), Gu = Ainv),
    rcov   = ~ vsm(usm(trait), ism(record)),
    data = long_sommer,
    control = mmes_control(diagnostics = list(verbose = FALSE))))
# covmatrix_mmes() returns the trait covariance of each term: 1 is the animal
# effect (G) and 2 the residual (E)
results$sommer <- compare("sommer",
                          covmatrix_mmes(fit_sommer, 1)$covariance,
                          covmatrix_mmes(fit_sommer, 2)$covariance,
                          elapsed = time_sommer[["elapsed"]])

## lme4breeding: REML through lme4, relationship matrix in relmat() ----------

if (requireNamespace("lme4breeding", quietly = TRUE)) {
  # lmeb() calls lme4 helpers from the search path, so attach the package
  library(lme4breeding)
  long <- data.frame(
    id = factor(rep(dat$id, n_traits)),
    trait = factor(rep(traits, each = nrow(dat)), levels = traits),
    value = unlist(dat[traits], use.names = FALSE),
    sex = rep(dat$sex, n_traits),
    age = rep(dat$age, n_traits))
  long$units <- long$id
  time_l4b <- system.time(
    fit_l4b <- lmeb(
      value ~ 0 + trait + trait:sex + trait:age +
        (0 + trait | id) + (0 + trait | units),
      relmat = list(id = A), data = long,
      REML = TRUE, rotation = FALSE, verbose = 0, trace = 0))
  vc_l4b <- VarCorr(fit_l4b)
  E_l4b <- as.matrix(vc_l4b$units)
  diag(E_l4b) <- diag(E_l4b) + attr(vc_l4b, "sc")^2
  results$lme4breeding <- compare("lme4breeding", as.matrix(vc_l4b$id),
                                  E_l4b, elapsed = time_l4b[["elapsed"]])
} else {
  message("lme4breeding is not installed; install.packages(\"lme4breeding\") ",
          "to add it")
}

## ASReml-R: REML, inverse of A passed through vm() --------------------------

if (!has_asreml) {
  message("asreml is not installed; it is licensed and not on CRAN, so the ",
          "ASReml-R fit is skipped")
} else {
  library(asreml)
  # Start from identity covariances, in the row-wise lower-triangle order us()
  # expects. The first pass can report non-convergence; the update finishes it.
  init_cov <- function(k) {
    m <- diag(k)
    m[cbind(rep(seq_len(k), seq_len(k)), sequence(seq_len(k)))]
  }
  Ainv_asr <- sparse_ginv(A)
  time_asr <- system.time({
    fit_asr <- suppressWarnings(asreml(
      fixed  = as.formula(paste(response, "~ trait + trait:sex + trait:age")),
      random = ~ vm(id, Ainv_asr):us(trait, init = init_cov(n_traits)),
      rcov   = ~ units:us(trait, init = init_cov(n_traits)),
      data = dat))
    fit_asr <- update(fit_asr)
  })
  vc <- summary(fit_asr)$varcomp
  results$asreml <- compare(
    "asreml",
    lower2mat(vc[grep("^vm\\(", rownames(vc)), 1], n_traits),
    lower2mat(vc[grep("^units.*trait!trait", rownames(vc)), 1], n_traits),
    elapsed = time_asr[["elapsed"]])
}

## The comparison -----------------------------------------------------------

results <- do.call(rbind, results)
rownames(results) <- results$engine
out <- round(results[, -1], 4)
out[, "time_s"] <- round(results$time_s, 1)
cat("\n", n_traits, "traits ,", n, "individuals \n", sep = "")
print(out)

# The four relative MSE columns divide each elementwise mean squared error by
# the mean squared truth on the same scale, so they are comparable both across
# quantities and across the number of traits; corrG uses the upper triangle,
# since the diagonal is always 1. min_ess is the smallest
# bulk effective sample size across the unique elements of G, and ess_ratio
# divides it by the retained draws (4000 for stanimal, 1000 for MCMCglmm and
# BGLR), so it measures mixing rather than chain length. REML engines have no
# draws and leave both columns NA. Timings are a single run each, in seconds,
# and cover the fit call only. Setup outside the call is not timed:
# animal_design() builds A for stanimal, solve(A) inverts it for sommer, and
# sparse_ginv(A) builds ASReml-R's inverse before its fit, while the inverse of A
# that MCMCglmm needs is evaluated inside its call and so is included.
# stanimal's time also includes compiling the model when the cache is cold;
# after the first fit the compiled executable is reused and only sampling is
# timed. These are not like-for-like posteriors either: MCMCglmm, BGLR and
# stanimal sample and stop at different effective sample sizes, while sommer,
# lme4breeding and ASReml-R maximize the REML likelihood.
