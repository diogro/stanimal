# Calomys skull and mandible data as fitted by stanimal
#
# The raw objects live in
# Dropbox/labbio/data/Calomys/{skull,skull2,mand}.RData and were assembled for a
# cross-classified sire and dam analysis (see main.R in that folder). Each trait
# is measured twice and the two replicas are averaged here.
#
# The skull set is the same 35 landmark distances used in the mice, so the module
# classification of Porto et al. (2009) applies without change; the trait names
# are converted from ISPM to IS_PM so they match ratones_morpho.
#
# The design is a half-sib design: 21 sires each mated to one to five of 16 dams,
# giving 55 full-sib cells. Parents are not themselves measured, so they are added
# to the pedigree as founders, which is what makeA() needs.
#
# Run from the repository root:
#   Rscript data-raw/calomys.R
# Writes data-raw/calomys.rds.

raw_dir <- "/Users/diogro/Library/CloudStorage/Dropbox/labbio/data/Calomys"

load(file.path(raw_dir, "skull.RData"))   # exp, exp2, main
load(file.path(raw_dir, "skull2.RData"))  # exp2
load(file.path(raw_dir, "mand.RData"))    # main

underscore <- function(x) {
  known <- c("IS", "NSL", "PNS", "PM", "ZS", "ZI", "MT", "NA", "BR", "PT",
             "APET", "BA", "EAM", "ZYGO", "TSP", "LD", "AS", "OPI", "JP")
  vapply(x, function(trait) {
    for (left in known)
      if (startsWith(trait, left) && nchar(trait) > nchar(left))
        return(paste0(left, "_", substring(trait, nchar(left) + 1L)))
    trait
  }, character(1), USE.NAMES = FALSE)
}

# The metadata columns are not in the same order in the two experiments, so they
# are taken by name rather than by position.
info_of <- function(x) {
  out <- x[, c("num", "sex", "age", "sir", "dam", "lit", "sib"), drop = FALSE]
  out$num <- as.character(out$num)
  out$sir <- as.character(out$sir)
  out$dam <- as.character(out$dam)
  as.data.frame(out, stringsAsFactors = FALSE)
}

average_replicas <- function(a, b) {
  traits <- colnames(a)[8:42]
  out <- (a[, traits] + b[, traits]) / 2
  colnames(out) <- underscore(traits)
  out
}

skull_traits <- rbind(average_replicas(exp[[1]], exp[[2]]),
                      average_replicas(exp2[[1]], exp2[[2]]))
skull_info <- rbind(info_of(exp[[1]]), info_of(exp2[[1]]))
skull_info$experiment <- rep(c("first", "second"), c(nrow(exp[[1]]), nrow(exp2[[1]])))

# Two individuals were measured in both experiments. Keeping the first record
# avoids a duplicated id, which makeA() would refuse.
keep <- !duplicated(skull_info$num)

skull <- data.frame(
  id = skull_info$num[keep],
  sire = skull_info$sir[keep],
  dam = skull_info$dam[keep],
  sex = factor(skull_info$sex[keep]),
  age = as.numeric(skull_info$age[keep]),
  litter = factor(skull_info$lit[keep]),
  experiment = factor(skull_info$experiment[keep]),
  stringsAsFactors = FALSE)
skull[names(skull_traits)] <- log(as.matrix(skull_traits[keep, ]))

# Mandible: columns 8 to 32 are the 25 distances and column 33 ("me") is a
# family code, not a measurement. main.R then drops every distance through the
# ramus tip, which leaves 20.
mand_traits <- main[, 8:32]
mand_traits <- mand_traits[, -grep("RT", colnames(mand_traits))]
mand <- data.frame(
  id = as.character(main$animal),
  sire = as.character(main$sir),
  dam = as.character(main$dam),
  sex = factor(main$sex),
  age = as.numeric(main$age),
  litter = factor(main$lit),
  family = factor(main$me),
  stringsAsFactors = FALSE)
mand[colnames(mand_traits)] <- log(as.matrix(mand_traits))

# Pedigree: the measured individuals plus every parent as a founder, since the
# parents are not themselves in the data.
pedigree_of <- function(d) {
  parents <- unique(c(d$sire, d$dam))
  rbind(
    data.frame(id = d$id, sire = d$sire, dam = d$dam, stringsAsFactors = FALSE),
    data.frame(id = setdiff(parents, d$id), sire = NA_character_,
               dam = NA_character_, stringsAsFactors = FALSE))
}

skull_pedigree <- pedigree_of(skull)
mand_pedigree <- pedigree_of(mand)

skull_trait_names <- names(skull)[grepl("_", names(skull))]
mand_trait_names <- setdiff(names(mand), c("id", "sire", "dam", "sex", "age",
                                          "litter", "family"))
cat("skull:", nrow(skull), "individuals,", length(skull_trait_names), "traits,",
    length(unique(skull$sire)), "sires,", length(unique(skull$dam)), "dams,",
    nrow(skull_pedigree), "pedigree rows\n")
cat("mandible:", nrow(mand), "individuals,", length(mand_trait_names),
    "traits,", length(unique(mand$sire)), "sires,", length(unique(mand$dam)),
    "dams,", nrow(mand_pedigree), "pedigree rows\n")

saveRDS(list(skull = skull, skull_pedigree = skull_pedigree,
             mandible = mand, mandible_pedigree = mand_pedigree),
        "data-raw/calomys.rds")
cat("Wrote data-raw/calomys.rds\n")
