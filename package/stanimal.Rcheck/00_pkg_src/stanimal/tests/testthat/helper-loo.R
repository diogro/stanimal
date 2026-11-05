skip_if_no_loo <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("loo")
}

# A small fit is enough for the bookkeeping; the checks are about shapes and
# plumbing, not about the posterior.
small_animal_fit <- function(...) {
  sim <- simulate_animal(n_fam = 8, size = 5, G = G3, E = E3, seed = 71)
  fit <- suppressWarnings(lmm_animal(sim$Y, sim$X, sim$A, iter = 400,
                                     chains = 2, cores = 2, seed = 1,
                                     refresh = 0, show_messages = FALSE, ...))
  list(fit = fit, sim = sim)
}
