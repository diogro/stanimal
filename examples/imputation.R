# Incomplete records in an animal model
#
# Traits are rarely measured on everyone. Dropping the individuals with a gap
# throws away the traits that *were* measured on them, and with relatedness it is
# unnecessary: a missing value is informed by every relative's observed traits.
# The opossum data ship that way -- 22 animals have no trait measured at all --
# so the model is fitted with the gaps already in place and imputed_data() fills
# them from the relatives' records.
#
# Run from the project root:
#   Rscript examples/imputation.R

library(stanimal)
library(cmdstanr)   # install.packages("cmdstanr",
                    #   repos = c("https://stan-dev.r-universe.dev", getOption("repositories")))
library(ggplot2)
options(mc.cores = 4)

## Data, with the gaps already in it
data(opossum)

# Four distances spanning the skull from the nasal suture to the basicranium,
# positively correlated with each other (0.39 to 0.78), which is what lets a gap
# borrow from the other traits of the same animal on top of what its relatives
# say. 11 of the animals whose sex and age are recorded have none of the four.
traits <- c("IS_NSL", "NA_BR", "ZI_ZYGO", "BA_EAM")
dat <- opossum
dat[traits] <- log(dat[traits])   # morphometrics are fitted on the log scale

# animal_design() drops the rows with a missing predictor (age), but keeps the
# rows whose response has gaps -- that is the whole point here.
design <- suppressWarnings(
  animal_design(cbind(IS_NSL, NA_BR, ZI_ZYGO, BA_EAM) ~ sex + age, data = dat))
Y <- design$Y
A <- opossum_grm[rownames(Y), rownames(Y)]
cat(sum(is.na(Y)), "of", length(Y), "entries are missing, on",
    sum(!complete.cases(Y)), "animals\n")
cat("dropping incomplete rows would leave", sum(complete.cases(Y)), "of",
    nrow(Y), "animals\n\n")

## Fit with the gaps in place: nothing is dropped
fit <- lmm_animal(Y, design$X, A, lkj_prior = 2,
                  iter = 800, warmup = 400, chains = 2)

## The completed data
imp <- imputed_data(fit, ndraws = 500, seed = 2)
imp                                # each gap: posterior mean and 95% interval

## Plot: what the model believes about each gap
gaps <- imp$gaps
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
  hi = intervals[, "hi"]
)
plot_dat$gap <- factor(plot_dat$gap, levels = rev(plot_dat$gap[order(plot_dat$trait)]))

p <- ggplot(plot_dat, aes(x = mean, y = gap)) +
  geom_pointrange(aes(xmin = lo, xmax = hi), colour = "grey30") +
  facet_wrap(~trait, scales = "free", nrow = 1) +
  labs(x = "log distance", y = NULL,
       title = "Imputed values for the animals with no measurement",
       subtitle = "posterior mean and 95% interval for each missing value") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 7))

dir.create("examples/figures", showWarnings = FALSE)
ggsave("examples/figures/imputation.png", p, width = 10, height = 5, dpi = 150)
cat("\nwrote examples/figures/imputation.png\n")
