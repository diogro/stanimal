ped6 <- data.frame(id = 1:6, sire = c(0, 0, 1, 1, 4, 4),
                   dam = c(0, 0, 2, 2, 3, 5))

# The tabular method, worked out by hand for ped6: individuals 5 and 6 are
# inbred through the full-sib mating of 4 with 3 and of 5 with 4.
expected_ped6 <- structure(
  c(1, 0, 0.5, 0.5, 0.5, 0.5,
    0, 1, 0.5, 0.5, 0.5, 0.5,
    0.5, 0.5, 1, 0.5, 0.75, 0.625,
    0.5, 0.5, 0.5, 1, 0.75, 0.875,
    0.5, 0.5, 0.75, 0.75, 1.25, 1,
    0.5, 0.5, 0.625, 0.875, 1, 1.375),
  dim = c(6L, 6L),
  dimnames = list(as.character(1:6), as.character(1:6)))

test_that("makeA computes the tabular relationship matrix", {
  A <- makeA(ped6)
  expect_equal(A, expected_ped6)
  expect_true(isSymmetric(A))
  expect_equal(unname(diag(A)), c(1, 1, 1, 1, 1.25, 1.375))
  expect_gt(min(eigen(A, symmetric = TRUE, only.values = TRUE)$values), 0)
})

test_that("the order of the pedigree does not matter", {
  shuffled <- ped6[c(5, 1, 6, 3, 2, 4), ]
  A <- makeA(shuffled)
  expect_equal(unname(A[rownames(expected_ped6), colnames(expected_ped6)]),
               unname(expected_ped6))
})

test_that("unknown parents can be coded 0, NA or empty", {
  from_na <- data.frame(id = 1:6, sire = c(NA, NA, 1, 1, 4, 4),
                        dam = c(NA, NA, 2, 2, 3, 5))
  from_empty <- ped6
  from_empty$sire <- as.character(from_empty$sire)
  from_empty$dam <- as.character(from_empty$dam)
  from_empty$sire[1:2] <- ""
  from_empty$dam[1:2] <- ""
  expect_equal(makeA(from_na), expected_ped6)
  expect_equal(makeA(from_empty), expected_ped6)
})

test_that("character IDs are kept", {
  ped <- data.frame(id = letters[1:4], sire = c(NA, NA, "a", "c"),
                    dam = c(NA, NA, "b", "b"))
  A <- makeA(ped)
  expect_identical(rownames(A), letters[1:4])
  expect_identical(colnames(A), letters[1:4])
  expect_equal(unname(A["c", "d"]), 0.75)
  expect_equal(unname(A["a", "d"]), 0.25)
})

test_that("a parent that is not in the pedigree is added as a founder", {
  ped <- data.frame(id = 1:3, sire = c(NA, 1, 9), dam = c(NA, NA, NA))
  expect_warning(A <- makeA(ped), "added as founders: 9")
  expect_identical(rownames(A), c("1", "2", "3", "9"))
  expect_equal(unname(A["3", "9"]), 0.5)
  expect_equal(unname(A["3", "1"]), 0)
})

test_that("makeA accepts related founders through Gamma", {
  Gamma <- matrix(c(1, 0.5, 0.5, 1), 2,
                  dimnames = list(c("1", "2"), c("1", "2")))
  A <- makeA(ped6, Gamma = Gamma)
  expect_equal(unname(A["1", "2"]), 0.5)
  expect_equal(unname(A["3", "3"]), 1.25)   # 1 + 0.5 * Gamma[1, 2]
  expect_equal(unname(A["3", "1"]), 0.75)
  expect_equal(unname(A["5", "5"]), 1.375)
  expect_true(isSymmetric(A))
  expect_gt(min(eigen(A, symmetric = TRUE, only.values = TRUE)$values), 0)

  identity <- diag(2)
  dimnames(identity) <- list(c("1", "2"), c("1", "2"))
  expect_equal(makeA(ped6, Gamma = identity), makeA(ped6))

  # the founder diagonal can differ from 1 in a metafounder model
  Gamma2 <- matrix(c(1.5, 0.5, 0.5, 1.5), 2,
                   dimnames = list(c("1", "2"), c("1", "2")))
  A2 <- makeA(ped6, Gamma = Gamma2)
  expect_equal(unname(diag(A2)[1:2]), c(1.5, 1.5))
  expect_equal(unname(A2["3", "3"]), 1.25)
  expect_equal(unname(A2["3", "1"]), 1)
})

test_that("Gamma also covers parents added as founders", {
  ped <- data.frame(id = 1:3, sire = c(NA, 1, 9), dam = c(NA, NA, NA))
  Gamma <- matrix(c(1, 0.25, 0.25, 1), 2,
                  dimnames = list(c("1", "9"), c("1", "9")))
  expect_warning(A <- makeA(ped, Gamma = Gamma), "added as founders: 9")
  expect_equal(unname(A["1", "9"]), 0.25)
  expect_equal(unname(A["2", "9"]), 0.125)
  expect_equal(unname(A["3", "9"]), 0.5)
})

test_that("makeA accepts metafounder groups", {
  Gamma <- matrix(c(1, 0.2, 0.2, 1), 2,
                  dimnames = list(c("A", "B"), c("A", "B")))
  ped <- data.frame(id = 1:5, sire = c(NA, NA, NA, 1, 1),
                    dam = c(NA, NA, NA, 2, 3))
  A <- makeA(ped, Gamma = Gamma, groups = c("A", "A", "B"))
  expect_equal(unname(A["1", "2"]), 1)      # same group
  expect_equal(unname(A["1", "3"]), 0.2)    # different groups
  expect_equal(unname(A["4", "4"]), 1.5)    # 1 + 0.5 * Gamma[A, A]
  expect_equal(unname(A["5", "5"]), 1.1)    # 1 + 0.5 * Gamma[A, B]

  # group and founder order can be given by names
  named_groups <- c(`2` = "A", `1` = "A", `3` = "B")
  expect_equal(makeA(ped, Gamma = Gamma, groups = named_groups), A)
  expect_equal(makeA(ped, Gamma = unname(Gamma), groups = c("A", "A", "B")), A)
})

test_that("makeA validates Gamma and groups", {
  expect_error(makeA(ped6, groups = c("A", "A")), "requires `Gamma`")
  expect_error(makeA(ped6, Gamma = c(1, 0, 0, 1)), "finite numeric matrix")
  expect_error(makeA(ped6, Gamma = diag(3)), "one row and column per founder")
  expect_error(makeA(ped6, Gamma = matrix(c(1, 2, 0, 1), 2)), "symmetric")
  expect_error(makeA(ped6, Gamma = matrix(c(1, 2, 2, 1), 2)), "positive definite")
  expect_error(makeA(ped6, Gamma = diag(2), groups = "A"),
               "one entry per founder")
  expect_error(makeA(ped6, Gamma = diag(2), groups = c("A", NA)),
               "missing values")
  Gamma_named <- matrix(c(1, 0, 0, 1), 2,
                        dimnames = list(c("A", "B"), c("A", "B")))
  expect_error(makeA(ped6, Gamma = Gamma_named, groups = c("C", "D")),
               "no row for group")
  Gamma_missing <- diag(2)
  dimnames(Gamma_missing) <- list(c("1", "9"), c("1", "9"))
  expect_error(makeA(ped6, Gamma = Gamma_missing), "no row for founder")
})

# nadiv is not needed to use the package: it is suggested only so that
# makeA() can be checked against an independent implementation.
test_that("makeA matches nadiv::makeA", {
  skip_if_not_installed("nadiv")
  compare <- function(pedigree) {
    reference <- as.matrix(suppressWarnings(nadiv::makeA(nadiv::prepPed(pedigree))))
    got <- suppressWarnings(makeA(pedigree))
    expect_equal(got, reference[rownames(got), colnames(got)], tolerance = 1e-12)
  }

  compare(ped6)                                             # the test pedigree
  compare(ped6[c(5, 1, 6, 3, 2, 4), ])                      # unordered
  compare(data.frame(id = 1:6, sire = c(NA, NA, NA, 1, 1, 4),  # full-sib mating
                     dam = c(NA, NA, NA, 2, 2, 5)))
  compare(data.frame(id = 3:1, sire = c(NA, 3, NA),            # parent added as founder
                     dam = c(NA, NA, NA)))
  compare(data.frame(id = letters[1:4], sire = c(NA, NA, "a", "c"),
                     dam = c(NA, NA, "b", "b")))
  compare(ratones_pedigree)                                 # 525 individuals

  set.seed(4)
  n <- 200
  founders <- 20
  random <- data.frame(id = seq_len(n),
                       sire = c(rep(NA, founders),
                                sample.int(founders, n - founders, TRUE)),
                       dam = c(rep(NA, founders),
                               sample.int(founders, n - founders, TRUE)))
  compare(random)
})
