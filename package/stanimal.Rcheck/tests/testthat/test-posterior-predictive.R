# The tail probability is the share of replicates at or above the observation,
# which is a property of the definition and can be checked without a sampler: an
# observation sitting at the centre of the replicates must land near 0.5, and one
# outside them must land at the ends.
test_that("the predictive tail probability counts replicates above the data", {
  expect_equal(stanimal:::.predictive_tail(c(1, 2, 3, 4, 5), 3), 0.6)
  expect_equal(stanimal:::.predictive_tail(c(1, 2, 3), 10), 0)
  expect_equal(stanimal:::.predictive_tail(c(1, 2, 3), 0), 1)

  set.seed(1)
  replicated <- stats::rnorm(1000)
  expect_equal(stanimal:::.predictive_tail(replicated, stats::median(replicated)),
               0.5, tolerance = 0.05)
  # a feature the model cannot produce is at the ends, not in the middle
  expect_lt(stanimal:::.predictive_tail(replicated,
                                        stats::quantile(replicated, 0.999)), 0.01)
})

test_that("a well-specified fit does not flag the features it reproduces", {
  skip_if_no_cmdstan()
  sim <- simulate_animal(n_fam = 12, size = 7, G = G3, E = E3, seed = 81)
  fit <- suppressWarnings(lmm_animal(sim$Y, sim$X, sim$A, iter = 600,
                                     chains = 2, cores = 2, seed = 1,
                                     refresh = 0, show_messages = FALSE))
  pp <- posterior_predictive(fit, ndraws = 100, seed = 2)
  out <- capture.output(print(pp))
  # the table's p column: parse it rather than re-deriving it, so the test
  # covers what the user reads
  rows <- grep("^ +(mean|sd|cor): ", out, value = TRUE)
  expect_length(rows, 3L + 3L + 3L)
  p <- as.numeric(sub("^.* ([0-9.]+)$", "\\1", rows))
  # under a correct model these are roughly uniform, so no single value is
  # informative: what a good fit cannot produce is p at the very ends
  expect_equal(sum(p == 0 | p == 1), 0L)
  expect_gte(sum(p > 0.1 & p < 0.9), 5L)
})
