# Breeding values from an animal model
#
# An animal model answers "how much of the variation in this trait is additive"
# and gives the covariance G. It does not answer "what is the additive value of
# *this* mouse", because the breeding values are integrated out of the
# likelihood: that is what makes the model scale to many individuals and many
# traits. breeding_values() puts them back afterwards, from the G, E and beta the
# fit already estimated. What comes out is what a breeding program would act on:
# a predicted additive value per individual, with the reliability of that
# prediction.
#
# Run from the project root:
#   Rscript examples/breeding_values.R

library(stanimal)
library(cmdstanr)   # install.packages("cmdstanr",
                    #   repos = c("https://stan-dev.r-universe.dev", getOption("repositories")))
library(ggplot2)
options(mc.cores = 4)

## Data
data(ratones_morpho)
# Three distances that span the skull and are not dominated by a single axis.
traits <- c("IS_NSL", "NA_BR", "ZI_ZYGO")
dat <- ratones_morpho
dat[traits] <- log(dat[traits])   # morphometrics are fitted on the log scale

# Each mouse records its dam and sire, so the pedigree is rebuilt from the data
# itself; parents that were never measured are added as founders by makeA().
ped <- data.frame(id = as.character(dat$ID),
                  sire = as.character(dat$PADRE),
                  dam = as.character(dat$MADRE),
                  stringsAsFactors = FALSE)

design <- suppressWarnings(
  animal_design(cbind(IS_NSL, NA_BR, ZI_ZYGO) ~ SEX + AGE + LIN,
                dat, pedigree = ped))

fit <- lmm_animal(design$Y, design$X, design$A, lkj_prior = 2,
                  iter = 800, warmup = 400, chains = 2)

## The breeding values
ebv <- breeding_values(fit)
ebv                       # the first mice, with standard errors and reliabilities

## Reliability: the share of a trait's additive variance the prediction pins
## down. A mouse with many measured relatives scores higher than a lone one.
cat("\nreliability (5%, 50%, 95% of mice):\n")
print(round(apply(ebv$reliability, 2, stats::quantile, c(0.05, 0.5, 0.95)), 2))

## The prediction is on the trait's scale, as a deviation from the fixed effects
residual <- design$Y - design$X %*% t(draws_mean(fit, "beta"))
cat("\nrank correlation between the breeding value and the mouse's own",
    "deviation:\n")
print(round(diag(stats::cor(ebv$values, residual, method = "spearman")), 2))

## Ranking candidates: the two orders are not the same list
first <- traits[1]
by_value <- rownames(ebv$values)[order(ebv$values[, first], decreasing = TRUE)]
by_phenotype <- rownames(residual)[order(residual[, first], decreasing = TRUE)]
cat("\ntop 10 mice by breeding value, with reliability:\n")
print(data.frame(id = by_value[1:10],
                 value = round(ebv$values[by_value[1:10], first], 3),
                 se = round(ebv$se[by_value[1:10], first], 3),
                 reliability = round(ebv$reliability[by_value[1:10], first], 3)),
      row.names = FALSE)
cat("\n", length(intersect(by_value[1:10], by_phenotype[1:10])), "of those",
    "are also in the top 10 by their own measurement; the rest are promoted by",
    "their relatives\n")

## Shrinkage: every prediction sits closer to the mean than the mouse's own
## deviation, and the points that sit furthest inside the line are the ones with
## the lowest reliability
plot_dat <- data.frame(
  trait = rep(traits, each = nrow(residual)),
  measured = as.numeric(residual),
  predicted = as.numeric(ebv$values),
  reliability = as.numeric(ebv$reliability)
)
# log distances have their own scale per trait, so both axes are on the sd scale
plot_dat$measured <- plot_dat$measured / rep(sqrt(diag(draws_mean(fit, "P"))),
                                             each = nrow(residual))
plot_dat$predicted <- plot_dat$predicted / rep(sqrt(diag(draws_mean(fit, "P"))),
                                               each = nrow(residual))
plot_dat$trait <- factor(plot_dat$trait, levels = traits)

p <- ggplot(plot_dat, aes(x = measured, y = predicted, colour = reliability)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 1.4, alpha = 0.8) +
  scale_colour_viridis_c(option = "magma", direction = -1) +
  facet_wrap(~trait, scales = "free") +
  labs(x = "measured deviation from the fixed effects (phenotypic sds)",
       y = "predicted breeding value (phenotypic sds)",
       colour = "reliability",
       title = "Breeding values shrink the measurements toward the mean",
       subtitle = paste("dashed line: the measurement itself; the further inside",
                        "the line, the less the mouse's own record is trusted")) +
  theme_minimal()

dir.create("examples/figures", showWarnings = FALSE)
ggsave("examples/figures/breeding_values.png", p, width = 9, height = 4.2,
       dpi = 150)
cat("\nwrote examples/figures/breeding_values.png\n")
