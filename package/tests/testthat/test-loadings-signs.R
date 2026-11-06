# Loadings are identified only up to sign, so the shaped draws are aligned before
# they are summarised. These tests use hand-built draws, where the right answer is
# known, so they do not depend on a sampler.

flipping_loadings <- function() {
  # two traits, one factor, ten draws: the factor flips sign half way through, so
  # the unaligned posterior mean of every loading is zero. h2 is there because
  # fit_draws() reads the trait count from the animal-model components.
  signs <- c(rep(1, 5), rep(-1, 5))
  cbind("h2[1]" = 0.5, "h2[2]" = 0.5,
        "loadings[1,1]" = 0.8 * signs,
        "loadings[2,1]" = 0.4 * signs)
}

test_that("sign switching no longer cancels a factor out of the mean", {
  m <- flipping_loadings()
  d <- fit_draws(m, traits = c("t1", "t2"))

  expect_equal(dim(d$loadings), c(2L, 1L, 10L))
  # every draw now has the same orientation, so the mean is the loading and not 0
  expect_true(all(d$loadings[, 1, ] > 0))
  expect_equal(unname(draws_mean(d, "loadings")[, 1]), c(0.8, 0.4))
  expect_equal(unname(draws_quantile(d, "loadings")[1, 1, ]), c(0.8, 0.8))
})

test_that("relabeling factor draws changes trait names but not factor names", {
  for (m in list(flipping_loadings(),
                 cbind(flipping_loadings(), "loadings[1,2]" = 0.2,
                       "loadings[2,2]" = 0.6,
                       "factor_magnitude[1]" = 1, "factor_magnitude[2]" = 2))) {
    d <- fit_draws(m, traits = c("old1", "old2"))
    original <- d
    means <- draws_mean(d, traits = c("new1", "new2"))
    q <- draws_quantile(d, "loadings", traits = c("new1", "new2"))
    expect_equal(rownames(means$loadings), c("new1", "new2"))
    expect_equal(dimnames(q)[[1]], c("new1", "new2"))
    expect_equal(colnames(means$loadings), dimnames(d$loadings)[[2]])
    expect_equal(dimnames(q)[[2]], dimnames(d$loadings)[[2]])
    expect_equal(unname(means$loadings), unname(draws_mean(d, "loadings")))
    expect_equal(unname(q), unname(draws_quantile(d, "loadings")))
    expect_equal(names(means$h2), c("new1", "new2"))
    if (!is.null(d$factor_magnitude))
      expect_equal(names(means$factor_magnitude), rownames(d$factor_magnitude))
    expect_identical(d, original)
  }
})

test_that("the anchor is chosen from the whole posterior, not per draw", {
  # factor 1 loads 0.9 on t1 and 0.95 on t3, and the two swap rank between
  # draws; a per-draw anchor would flip them differently
  m <- cbind("h2[1]" = 0.5, "h2[2]" = 0.5, "h2[3]" = 0.5,
             "loadings[1,1]" = c(0.9, 0.9, 0.9),
             "loadings[3,1]" = c(0.95, -0.95, 0.95),
             "loadings[2,1]" = c(0.5, -0.5, 0.5))
  d <- fit_draws(m, traits = c("t1", "t2", "t3"))
  # the largest mean absolute loading is t3, so t3 is positive in every draw and
  # draw 2 is the one that had to be flipped
  expect_true(all(d$loadings[3, 1, ] > 0))
  expect_equal(unname(d$loadings[1, 1, ]), c(0.9, -0.9, 0.9))
})

test_that("a factor that never switches is left alone", {
  m <- cbind("h2[1]" = 0.5, "h2[2]" = 0.5,
             "loadings[1,1]" = c(0.1, 0.2), "loadings[2,1]" = c(-0.3, -0.4),
             "loadings[2,2]" = c(0.5, 0.6))
  d <- fit_draws(m, traits = c("t1", "t2"))
  # factor 1 is anchored on t2 (|t2| is larger), so the whole column flips and
  # the relative pattern within it is preserved
  expect_equal(unname(d$loadings[1, 1, ]), c(-0.1, -0.2))
  expect_equal(unname(d$loadings[2, 1, ]), c(0.3, 0.4))
  # factor 2 was already positive, so nothing moves
  expect_equal(unname(d$loadings[2, 2, ]), c(0.5, 0.6))
})

test_that("the raw draws are still available unaligned", {
  skip_if_no_cmdstan()
  sim <- simulate_animal(n_fam = 8, size = 5, G = G3, E = E3, seed = 81)
  fit <- suppressWarnings(lmm_animal(sim$Y, sim$X, sim$A,
                                     genetic_structure = "factor",
                                     n_factors = 2, column_decay = 0.5,
                                     iter = 400, chains = 2, cores = 2, seed = 1,
                                     refresh = 0, show_messages = FALSE))
  aligned <- fit_draws(fit)$loadings
  raw <- as.matrix(fit$draws("loadings", format = "draws_matrix"))
  expect_equal(dim(aligned), c(ncol(sim$Y), 2L, nrow(raw)))
  # alignment may flip a column, so compare absolute values: those are unchanged
  expect_equal(abs(unname(aligned[1, 1, ])),
               abs(as.numeric(raw[, "loadings[1,1]"])))
})
