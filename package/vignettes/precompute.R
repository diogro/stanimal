# Knit the vignette sources (*.Rmd.orig) into the precomputed vignettes (*.Rmd).
#
# The vignettes fit Stan models, which take minutes and need CmdStan, so R CMD
# check and pkgdown never run them: they render the .Rmd files written here, in
# which every output and figure is already in place. Rerun this after changing
# a vignette or anything its output depends on, with the current package
# installed (R CMD INSTALL package). From package/vignettes:
#
#   Rscript precompute.R                 # every vignette
#   Rscript precompute.R univariate gibbs
#
# Each vignette is knitted in its own R process, so none sees another's objects.
args <- commandArgs(trailingOnly = TRUE)
sources <- if (length(args)) paste0(args, ".Rmd.orig") else
  list.files(pattern = "\\.Rmd\\.orig$")
missing <- sources[!file.exists(sources)]
if (length(missing)) stop("no such vignette source: ", paste(missing, collapse = ", "))
for (src in sources) {
  out <- sub("\\.orig$", "", src)
  message("knitting ", src, " -> ", out)
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    c("-e", shQuote(sprintf("knitr::knit('%s', '%s', quiet = TRUE)",
                                            src, out))))
  if (status != 0) stop("knitting ", src, " failed")
}
