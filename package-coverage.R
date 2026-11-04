## Measure test coverage of the stanimal package for the manuscript.
## Writes a per-file and a per-package summary to coverage-report.txt and saves
## the coverage object to coverage.rds so the table can be rebuilt without
## running the suite again.
## Run from the repository root: Rscript package-coverage.R

suppressMessages(library(covr))

cov <- package_coverage(
  path  = "package",
  type  = "tests",
  clean = TRUE,
  quiet = FALSE
)
saveRDS(cov, "coverage.rds")

lines <- as.data.frame(tally_coverage(cov, by = "line"))
lines$covered <- lines$value > 0

by_file <- do.call(rbind, lapply(split(lines, lines$filename), function(d) {
  data.frame(filename = d$filename[1],
             covered  = sum(d$covered),
             total    = nrow(d))
}))
by_file$pct <- 100 * by_file$covered / by_file$total
by_file <- by_file[order(by_file$pct), ]

out <- file("coverage-report.txt", "w")
writeLines("Lines covered / instrumented, least covered first:", out)
for (i in seq_len(nrow(by_file))) {
  writeLines(sprintf("  %-26s %5d / %5d  (%5.1f%%)",
                     by_file$filename[i], by_file$covered[i],
                     by_file$total[i], by_file$pct[i]), out)
}
writeLines("", out)
writeLines(sprintf("Package total: %.1f%% of %d instrumented lines",
                   percent_coverage(cov), nrow(lines)), out)
close(out)

cat(readLines("coverage-report.txt"), sep = "\n")
