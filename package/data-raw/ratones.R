# Builds the bundled mouse dataset (data/ratones_morpho.rda) from the companion
# data package diogro/ratones.
#
# The original records come from Penna et al. (2017), archived at Dryad
# (doi:10.5061/dryad.5gr8r) and packaged in https://github.com/diogro/ratones,
# which is MIT licensed. Set RATONES_SOURCE to the data/ directory of a
# checkout of that repository, or leave it unset to use the default location.
#
# Run from the package directory: Rscript data-raw/ratones.R

source_dir <- Sys.getenv("RATONES_SOURCE",
                         unset = file.path("..", "..", "ratones", "data"))
load(file.path(source_dir, "ratonesdf.RData"))

ratones_morpho <- ratonesdf
ratones_morpho$ID <- as.integer(ratones_morpho$ID)
ratones_morpho$GER <- as.integer(ratones_morpho$GER)
ratones_morpho$MADRE <- as.integer(ratones_morpho$MADRE)
ratones_morpho$PADRE <- as.integer(ratones_morpho$PADRE)
# Row names carry the individual IDs, so animal_design() can match a pedigree
# or a relationship matrix to these rows without any extra preparation.
rownames(ratones_morpho) <- as.character(ratones_morpho$ID)

stopifnot(
  identical(dim(ratones_morpho), c(329L, 46L)),
  !anyDuplicated(ratones_morpho$ID)
)

dir.create("data", showWarnings = FALSE)
save(ratones_morpho, file = "data/ratones_morpho.rda", compress = "xz", version = 2)
