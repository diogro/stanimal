test_that("ratones_morpho has the documented shape", {
  expect_s3_class(ratones_morpho, "data.frame")
  expect_identical(dim(ratones_morpho), c(329L, 46L))
  expect_identical(rownames(ratones_morpho), as.character(ratones_morpho$ID))
  expect_identical(colnames(ratones_morpho)[1:11],
                   c(".id", "ID", "SEX", "GER", "LIN", "MADRE", "PADRE",
                     "AGE", "P49", "line", "selection"))
  expect_length(names(ratones_morpho)[12:46], 35L)
  expect_false(anyNA(ratones_morpho))
  expect_setequal(unique(ratones_morpho$selection),
                  c("control", "downwards", "upwards"))
  expect_true(all(ratones_morpho$SEX %in% c("F", "M")))
})

test_that("opossum has the documented shape", {
  expect_s3_class(opossum, "data.frame")
  expect_identical(dim(opossum), c(553L, 38L))
  expect_identical(rownames(opossum), as.character(opossum$id))
  expect_identical(names(opossum)[1:3], c("id", "sex", "age"))
  expect_identical(anyDuplicated(opossum$id), 0L)
  expect_false(anyNA(opossum$id))
  expect_false(anyNA(opossum$sex))
  expect_true(all(opossum$sex %in% c("1", "2")))
  # the same 35 landmarks as the mice
  expect_setequal(setdiff(names(opossum), c("id", "sex", "age")),
                  names(ratones_morpho)[12:46])
  # 22 animals have no trait measured, and 505 have every trait and covariate
  traits <- setdiff(names(opossum), c("id", "sex", "age"))
  expect_identical(sum(rowSums(is.na(opossum[traits])) == length(traits)), 22L)
  expect_identical(sum(stats::complete.cases(opossum)), 505L)
})

test_that("opossum_grm matches the animals and the complete subset is usable", {
  expect_true(isSymmetric(opossum_grm))
  expect_identical(rownames(opossum_grm), as.character(opossum$id))
  expect_identical(colnames(opossum_grm), rownames(opossum_grm))
  expect_equal(range(diag(opossum_grm)), c(0.2228368, 0.5942065),
               tolerance = 1e-4)
  # the full matrix is only positive semidefinite; the analysed subset is not
  ok <- stats::complete.cases(opossum)
  A <- opossum_grm[ok, ok]
  expect_gt(min(eigen(A, symmetric = TRUE, only.values = TRUE)$values), 0)
})

test_that("the morphometrics are moderately integrated", {
  traits <- names(ratones_morpho)[12:46]
  Y <- log(as.matrix(ratones_morpho[, traits]))
  pc1_share <- pc1_variance_share(Y)
  expect_gt(pc1_share, 0.3)
  expect_lt(pc1_share, 0.8)
  expect_gt(mean(cor(Y)[upper.tri(diag(ncol(Y)))]), 0.2)
})
