# Incomplete records in an animal model
#
# Traits are rarely measured on everyone. Dropping the individuals with a gap
# throws away the traits that *were* measured on them, and with relatedness it is
# unnecessary: a missing value is informed by every relative's observed traits.
# Here a few entries of the mouse data are blanked out, the model is fitted with
# the gaps in place, and imputed_data() puts back what was removed -- which can
# be scored, because the deleted values are known.
#
# Run from the project root:
#   Rscript examples/imputation.R

library(stanimal)
library(cmdstanr)   # install.packages("cmdstanr",
                    #   repos = c("https://stan-dev.r-universe.dev", getOption("repositories")))
library(ggplot2)
options(mc.cores = 4)

## Data, and gaps to fill
data(ratones_morpho)
data(ratones_pedigree)
# Four distances spanning the skull from the nasal suture to the basicranium.
# All four are positively correlated with each other (0.18 to 0.49), which is
# what lets a gap borrow from the other traits of the same mouse on top of what
# its relatives say.
traits <- c("IS_NSL", "NA_BR", "ZI_ZYGO", "BA_EAM")
dat <- ratones_morpho
dat[traits] <- log(dat[traits])   # morphometrics are fitted on the log scale

design <- animal_design(cbind(IS_NSL, NA_BR, ZI_ZYGO, BA_EAM) ~ SEX + AGE + LIN,
                        dat, pedigree = ratones_pedigree)
Y <- design$Y
truth <- Y                         # kept only to score the imputation below

set.seed(1)
scattered <- cbind(sample(nrow(Y), 8), rep(seq_len(4), each = 2))
Y[scattered] <- NA
Y[200:205, "ZI_ZYGO"] <- NA        # one trait missing on a whole set of mice
cat(sum(is.na(Y)), "of", length(Y), "entries removed,",
    sum(!complete.cases(Y)), "mice affected\n")
cat("dropping incomplete rows would leave", sum(complete.cases(Y)), "of",
    nrow(Y), "mice\n\n")

## Fit with the gaps in place: nothing is dropped
fit <- lmm_animal(Y, design$X, design$A, lkj_prior = 2,
                  iter = 800, warmup = 400, chains = 2)

## The completed data
imp <- imputed_data(fit, ndraws = 500, seed = 2)
imp                                # each gap: posterior mean and 95% interval

gaps <- imp$gaps
posterior_mean <- apply(imp$Y, c(1L, 2L), mean)[cbind(gaps$individual, gaps$trait)]
posterior_sd <- apply(imp$Y, c(1L, 2L), sd)[cbind(gaps$individual, gaps$trait)]
removed <- truth[cbind(gaps$individual, gaps$trait)]
cat("\nof the", nrow(gaps), "removed values,",
    round(100 * mean(abs(posterior_mean - removed) < 2 * posterior_sd)),
    "% are within two posterior standard deviations of what was there\n")

## Plot: what the model believes about each gap
intervals <- t(vapply(seq_len(nrow(gaps)), function(i) {
  v <- imp$Y[gaps$individual[i], gaps$trait[i], ]
  c(mean = mean(v),
    lo = unname(quantile(v, 0.025)),
    hi = unname(quantile(v, 0.975)))
}, c(mean = 0, lo = 0, hi = 0)))
plot_dat <- data.frame(
  gap = sprintf("%s / %s", gaps$individual_label, gaps$trait_label),
  trait = gaps$trait_label,
  mean = intervals[, "mean"],
  lo = intervals[, "lo"],
  hi = intervals[, "hi"],
  removed = truth[cbind(gaps$individual, gaps$trait)]
)
plot_dat$gap <- factor(plot_dat$gap, levels = rev(plot_dat$gap[order(plot_dat$trait)]))

p <- ggplot(plot_dat, aes(x = mean, y = gap)) +
  geom_pointrange(aes(xmin = lo, xmax = hi), colour = "grey30") +
  geom_point(aes(x = removed), shape = 17, colour = "firebrick", size = 2.5) +
  facet_wrap(~trait, scales = "free", nrow = 1) +
  labs(x = "log distance", y = NULL,
       title = "Imputed values, and the data that were removed",
       subtitle = "point and 95% interval: the imputation; triangle: the value deleted") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 7))

dir.create("examples/figures", showWarnings = FALSE)
ggsave("examples/figures/imputation.png", p, width = 10, height = 5, dpi = 150)
cat("\nwrote examples/figures/imputation.png\n")
