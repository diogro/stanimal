# Build the pkgdown site into the repository's docs/ (served by GitHub Pages).
#
# This file and pkgdown/extra.scss are shared, unchanged, by stanimal, indri,
# grnland and evomod2. Edit one copy and copy it to the others; per-package
# settings (accent colour, navbar groups, reference sections) live in
# _pkgdown.yml.
#
# The home page is the repository README, so the site and GitHub show the same
# text. When the package sits in a subdirectory of the repository (README one
# level up), the README is copied to pkgdown/index.md with its relative links
# fixed for the site: package figures stay local, every other repository path
# points at GitHub. pkgdown/index.md is generated; edit the README instead.
#
# Usage, from anywhere:
#   Rscript <package>/pkgdown/build_site.R          # full rebuild
#   Rscript <package>/pkgdown/build_site.R --lazy   # only changed pages

args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
if (length(script) != 1) stop("run this file with Rscript")
pkg <- normalizePath(file.path(dirname(script), ".."))
lazy <- "--lazy" %in% commandArgs(trailingOnly = TRUE)

git <- function(...) {
  system2("git", c("-C", shQuote(pkg), ...), stdout = TRUE)
}
root <- normalizePath(git("rev-parse", "--show-toplevel"))
index <- file.path(pkg, "pkgdown", "index.md")

if (root != pkg) {
  sub_dir <- substring(pkg, nchar(root) + 2)
  urls <- trimws(strsplit(read.dcf(file.path(pkg, "DESCRIPTION"),
                                   "URL")[1, 1], ",")[[1]])
  repo <- sub("/$", "", grep("^https://github.com/", urls, value = TRUE)[1])
  if (is.na(repo)) stop("DESCRIPTION URL needs the GitHub repository")
  branch <- tryCatch(
    sub("^origin/", "", git("rev-parse", "--abbrev-ref", "origin/HEAD")),
    warning = function(w) git("branch", "--show-current"))
  blob <- paste0(repo, "/blob/", branch, "/")

  relink <- function(target) {
    if (grepl("^([a-z]+:|#|/)", target)) return(target)
    prefix <- paste0(sub_dir, "/")
    if (startsWith(target, prefix)) {
      local <- substring(target, nchar(prefix) + 1)
      if (startsWith(local, "man/figures/")) return(local)
    }
    paste0(blob, target)
  }
  rewrite <- function(line, pattern) {
    m <- gregexpr(pattern, line, perl = TRUE)
    regmatches(line, m) <- lapply(regmatches(line, m), function(hits) {
      vapply(hits, function(hit) {
        target <- sub(pattern, "\\2", hit, perl = TRUE)
        sub(target, relink(target), hit, fixed = TRUE)
      }, character(1))
    })
    line
  }

  text <- readLines(file.path(root, "README.md"), warn = FALSE)
  # pkgdown shows man/figures/logo.* beside the title: drop the README's copy
  # (and a <p> wrapped around it on lines of its own).
  for (i in rev(grep("<img[^>]*man/figures/logo\\.", text))) {
    drop <- i
    if (i > 1 && i < length(text) && grepl("^\\s*<p[^>]*>\\s*$", text[i - 1]) &&
        grepl("^\\s*</p>\\s*$", text[i + 1])) drop <- (i - 1):(i + 1)
    text <- text[-drop]
  }
  fence <- cumsum(grepl("^\\s*```", text)) %% 2 == 1
  for (i in which(!fence)) {
    text[i] <- rewrite(text[i], "(\\]\\()([^)\\s]+)(\\))")
    text[i] <- rewrite(text[i], "((?:src|href)=\")([^\"]+)(\")")
  }
  writeLines(c("<!-- Generated from the repository README by",
               "     pkgdown/build_site.R: edit README.md, not this file. -->",
               "", text), index)
} else if (file.exists(index)) {
  stop(index, " would replace README.md as the home page; remove it")
}

pkgdown::build_site(pkg, lazy = lazy, preview = FALSE)
