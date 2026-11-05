# Opossum cranial data from Porto et al. (2016)
#
# Porto, Schmelter, VandeBerg, Marroig and Cheverud (2016), "Evolution of the
# genotype-to-phenotype map and the cost of pleiotropy in mammals", Genetics
# 204:1601-1612, doi:10.1534/genetics.116.189431. Data archived on figshare as
# "Data - Genetics - Porto et al 2016", https://figshare.com/articles/4055961.
#
# 553 measured Monodelphis domestica from the VandeBerg colony, 35 cranial
# Euclidean distances on the same landmark set as the mice, plus age and sex, a
# kinship matrix over the colony, and the published 35 x 35 G and P matrices.
# The published matrices are the point: they let a stanimal fit be checked
# against an independent analysis instead of against another of our own fits.
#
# Run from the repository root:
#   Rscript data-raw/opossum.R

src <- "/Users/diogro/Downloads/4055961"

cranial <- read.csv(file.path(src, "phenotype_covariates", "cranial_data.csv"),
                    stringsAsFactors = FALSE)
kin <- read.csv(file.path(src, "phenotype_covariates", "kinship_ped.csv"),
                stringsAsFactors = FALSE)
names(kin) <- c("ego1", "ego2", "kin")

# The published matrices, for the scale check and later for validation.
unzip(file.path(src, "matrices.zip"), exdir = tempdir())
G_pub <- as.matrix(readxl::read_excel(file.path(tempdir(), "matrices",
                                               "G_matrix.xlsx"), col_names = FALSE))
P_pub <- as.matrix(readxl::read_excel(file.path(tempdir(), "matrices",
                                               "P_matrix.xlsx"), col_names = FALSE))

landmarks <- c("IS", "NSL", "PNS", "PM", "ZS", "ZI", "MT", "NA", "BR", "PT",
               "APET", "BA", "EAM", "ZYGO", "TSP", "LD", "AS", "OPI", "JP")
landmarks <- landmarks[order(-nchar(landmarks))]

# The column names carry a trailing replica digit (PMZS1) that the mice names do
# not, so it is stripped before the landmark split.
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
traits <- porto_name(measured)
stopifnot(!anyDuplicated(traits))
cat("traits mapped:", length(traits), " of which unmatched:",
    sum(!grepl("_", traits)), "\n")

# Scale check: the published P is compared with the variance of the raw and of
# the logged traits, which says which scale the published analysis used.
raw_var <- apply(cranial[measured], 2, var)
log_var <- apply(log(as.matrix(cranial[measured])), 2, var)
cat("correlation of diag(P) with var(raw):", round(cor(diag(P_pub), raw_var), 3),
    " with var(log):", round(cor(diag(P_pub), log_var), 3), "\n")
cat("diag(P) range:", round(range(diag(P_pub)), 4),
    " var(log) range:", round(range(log_var), 4), "\n")
cat("published h2: mean", round(mean(diag(G_pub) / diag(P_pub)), 3),
    " range", round(range(diag(G_pub) / diag(P_pub)), 3), "\n")

# The kinship file is long form over every pair; rebuild the matrix.
ids <- sort(unique(c(kin$ego1, kin$ego2)))
A <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
A[cbind(kin$ego1, kin$ego2)] <- kin$kin
A[cbind(kin$ego2, kin$ego1)] <- kin$kin
cat("kinship:", length(ids), "animals, rank", qr(A)$rank,
    " symmetric:", isTRUE(all.equal(A, t(A))), "\n")

keep <- cranial$Animal %in% ids
cat("measured animals in the kinship matrix:", sum(keep), "of", nrow(cranial), "\n")

dat <- data.frame(
  id = cranial$Animal[keep],
  sex = factor(cranial$SexBI[keep]),
  age = as.numeric(cranial$AgeBI[keep]),
  stringsAsFactors = FALSE)
dat[traits] <- log(as.matrix(cranial[keep, measured]))

A <- A[dat$id, dat$id]
cat("h2 upper bound from the eigenvalues of A (should be < 1):",
    round(max(eigen(A, symmetric = TRUE, only.values = TRUE)$values), 3), "\n")

saveRDS(list(data = dat, A = A, traits = traits,
             G_published = G_pub, P_published = P_pub),
        "data-raw/opossum.rds")
cat("Wrote data-raw/opossum.rds\n")
