# Opossum cranial data and genomic relationship matrix
#
# Builds package/data/opossum.rda from the archive of Porto, Schmelter, VandeBerg,
# Marroig and Cheverud (2016), "Evolution of the genotype-to-phenotype map and the
# cost of pleiotropy in mammals", Genetics 204:1601-1612,
# doi:10.1534/genetics.116.189431, deposited on figshare as "Data - Genetics -
# Porto et al 2016", https://figshare.com/articles/dataset/4055961.
#
# 553 Monodelphis domestica from the VandeBerg colony, measured for the same 35
# cranial Euclidean distances as the mice, with age and sex. The relationship
# matrix is genomic, from the 3,696 coded markers, because the shipped pedigree
# kinship has every pair related at about 0.77 and cannot separate additive
# variance from a shared colony effect.
#
# The traits are stored as measured: Porto et al. analysed these distances
# untransformed, which their P matrix confirms (diagonal correlation 0.834 with
# the raw variances against 0.046 with the logged ones).
#
# Run from the package directory:
#   Rscript data-raw/opossum.R

article <- 4055961
cache <- file.path(tempdir(), "opossum-figshare")
dir.create(cache, recursive = TRUE, showWarnings = FALSE)

# The download links are resolved from the figshare API rather than hard-coded,
# since the file identifiers are not part of the article number. Everything is
# cached in the session temporary directory, so a rebuild costs nothing.
fetch <- function(name) {
  dest <- file.path(cache, name)
  if (!file.exists(dest)) {
    listing <- jsonlite::fromJSON(sprintf(
      "https://api.figshare.com/v2/articles/%d", article))$files
    url <- listing$download_url[listing$name == name]
    if (length(url) != 1L)
      stop("figshare article ", article, " does not offer ", name,
           call. = FALSE)
    utils::download.file(url, dest, mode = "wb", quiet = TRUE)
  }
  dest
}

unzip(fetch("phenotype_covariates.zip"), exdir = cache)
unzip(fetch("coded_marker_data.zip"), exdir = cache)

cranial <- read.csv(file.path(cache, "phenotype_covariates", "cranial_data.csv"),
                    stringsAsFactors = FALSE)

landmarks <- c("IS", "NSL", "PNS", "PM", "ZS", "ZI", "MT", "NA", "BR", "PT",
               "APET", "BA", "EAM", "ZYGO", "TSP", "LD", "AS", "OPI", "JP")
landmarks <- landmarks[order(-nchar(landmarks))]

# The measurement columns carry a trailing replica digit (PMZS1) that the mice
# names do not, so it is stripped before the landmark split.
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
stopifnot(!anyDuplicated(traits), all(grepl("_", traits)))

opossum <- data.frame(
  id = as.character(cranial$Animal),
  sex = factor(cranial$SexBI),
  age = as.numeric(cranial$AgeBI),
  stringsAsFactors = FALSE)
opossum[traits] <- as.matrix(cranial[measured])
rownames(opossum) <- opossum$id

# Markers: the coded files carry no animal IDs, so animals are matched to the
# cranial file by row order. Both have 553 rows; that equality is the only check
# available, and the package documentation says so.
marker_files <- list.files(file.path(cache, "coded_marker_data"),
                           pattern = "SAS_ADD", full.names = TRUE)
stopifnot(length(marker_files) == 8L)
markers <- do.call(cbind, lapply(marker_files, function(f)
  as.matrix(read.csv(f, header = FALSE))))
stopifnot(nrow(markers) == nrow(opossum), all(markers %in% c(-1, 0, 1)))

# The coding is additive, -1/0/1, so markers + 1 counts copies of one allele.
# The matrix is VanRaden's (2008) first estimator, Z Z' / (2 sum p (1 - p)) with
# Z centred by twice the sample allele frequency, which puts it on the scale of
# a pedigree relationship matrix (mean diagonal about 1) and is what
# genomic_relationship(markers + 1) returns: every marker passes its default
# filters (minor allele frequency at least 0.044, no missing calls). In
# stanimal 0.8.1 and earlier the denominator was the number of markers, which
# shrank the matrix by a factor of 2.86 (diagonal near 0.37), so the G fitted on
# it was inflated by about that factor and the heritabilities with it.
freq <- (colMeans(markers) + 1) / 2
stopifnot(!anyNA(markers), min(pmin(freq, 1 - freq)) > 0.01)
Z <- sweep(markers, 2, 2 * freq - 1, "-")
opossum_grm <- tcrossprod(Z) / (2 * sum(freq * (1 - freq)))
opossum_grm <- (opossum_grm + t(opossum_grm)) / 2
# The animal IDs are the link between the matrix and the phenotype frame; nothing
# in the marker files supplies them.
dimnames(opossum_grm) <- list(opossum$id, opossum$id)

cat("opossum:", nrow(opossum), "individuals,", length(traits), "traits,",
    "complete records:", sum(stats::complete.cases(opossum[traits])), "\n")
cat("opossum_grm:", nrow(opossum_grm), "x", ncol(opossum_grm), "from",
    ncol(markers), "markers\n")
cat("GRM diagonal range:", round(range(diag(opossum_grm)), 3),
    " mean:", round(mean(diag(opossum_grm)), 3),
    " off-diagonal sd:", round(stats::sd(opossum_grm[lower.tri(opossum_grm)]), 4), "\n")

save(opossum, opossum_grm, file = "data/opossum.rda", compress = "xz", version = 2)
cat("Wrote data/opossum.rda\n")
