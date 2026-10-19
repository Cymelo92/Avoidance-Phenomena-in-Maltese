# ============================================================
# Arabic offset model: prediction plots
# ============================================================
# This script produces two figures from the fitted Arabic
# offset negative binomial model:
#
# 1. A continuous prediction plot showing model-implied
#    relative frequency across the observed similarity range.
#
# 2. A contrast plot showing the estimated change in predicted
#    relative frequency between low similarity (z = -1) and
#    high similarity (z = 1).
#
# In both cases, expected_count is fixed at 1 so that predictions
# are interpretable on a relative-frequency scale rather than as
# absolute counts tied to a specific independence expectation.
# ============================================================

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(forcats)
library(scales)
library(readr)
library(brms)

# ------------------------------------------------------------
# 1. Load fitted model and prepare base data
# ------------------------------------------------------------
# The fitted brms model contains the original data used for fitting.
# We recover that data here and standardize factor ordering so that
# plots are consistent and interpretable.

fit_arabic_offset <- read_rds("results/arabic_offset.rds")

pair_counts2 <- fit_arabic_offset$data %>%
  mutate(
    position = factor(position, levels = c("12", "23", "13"))
  )

# ------------------------------------------------------------
# 2. Common relabeling helpers
# ------------------------------------------------------------
# These labels are used in both plots so that terminology stays
# consistent throughout the paper.

position_levels_pretty <- c(
  "12" = "C1-C2 (adjacent)",
  "23" = "C2-C3 (adjacent)",
  "13" = "C1-C3 (non-adjacent)"
)

# ------------------------------------------------------------
# 3. Plot 1: continuous predictions across similarity
# ------------------------------------------------------------
# We create a dense prediction grid over the observed similarity
# range, crossed with all positions.
#
# expected_count is set to 1 for all rows. Because the fitted model
# contains offset(log(expected_count)), this makes the fitted mean
# directly interpretable as model-implied relative frequency.

grid_offset <- expand.grid(
  similarity_z = seq(
    min(pair_counts2$similarity_z, na.rm = TRUE),
    max(pair_counts2$similarity_z, na.rm = TRUE),
    length.out = 200
  ),
  position = levels(pair_counts2$position),
  expected_count = 1
) %>%
  as_tibble()

# Obtain posterior fitted values on the response scale.
# Because expected_count = 1, these can be interpreted as
# predicted relative frequency.

epred_offset <- fitted(
  fit_arabic_offset,
  newdata = grid_offset,
  re_formula = NA,
  scale = "response",
  summary = TRUE
) %>%
  as.data.frame()

# Combine prediction grid with posterior summaries and create
# readable labels for plotting.

plot_ratio <- bind_cols(grid_offset, epred_offset) %>%
  mutate(
    rel_freq = Estimate / expected_count,
    rel_freq_low = Q2.5 / expected_count,
    rel_freq_high = Q97.5 / expected_count,
    position_label = recode(position, !!!position_levels_pretty)
  ) %>%
  mutate(
    position_label = factor(
      position_label,
      levels = unname(position_levels_pretty)
    )
  )

# Continuous prediction plot.
# Grey ribbons show uncertainty.

p_ratio <- ggplot(
  plot_ratio,
  aes(
    x = similarity_z,
    y = rel_freq,
    ymin = rel_freq_low,
    ymax = rel_freq_high
  )
) +
  geom_ribbon(alpha = 0.15, fill = "grey70", colour = NA) +
  geom_line() +
  facet_grid(. ~ position_label) +
  labs(
    x = "Phonetic similarity (z-score)",
    y = "Predicted relative frequency",
    title = "Predicted relative frequency by phonetic similarity",
    subtitle = "Offset model. Expected count fixed at 1"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title.position = "plot",
    strip.background = element_rect(fill = "grey85"),
    panel.grid.minor = element_blank()
  )

print(p_ratio)

ggsave(
  filename = "arabic_offset_ratio_conditional.png",
  plot = p_ratio,
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 4. Plot 2: contrast between low and high similarity
# ------------------------------------------------------------
# This plot summarizes the effect of similarity on the response
# scale by comparing posterior predictions at:
#   low similarity  = z = -1
#   high similarity = z = 1
#
# Again, expected_count is fixed at 1 so that the predictions are
# directly interpretable as relative frequency.

grid_low <- expand.grid(
  similarity_z = -1,
  position = levels(pair_counts2$position),
  expected_count = 1
) %>%
  as_tibble()

grid_high <- grid_low %>%
  mutate(similarity_z = 1)

# Get full posterior expected values for both grids.
# We use posterior_epred() rather than fitted() because we want
# draw-by-draw differences, not just posterior summaries.

epred_low <- posterior_epred(
  fit_arabic_offset,
  newdata = grid_low,
  re_formula = NA
)

epred_high <- posterior_epred(
  fit_arabic_offset,
  newdata = grid_high,
  re_formula = NA
)

# Compute posterior draws for the change in predicted relative
# frequency from low to high similarity.
#
# Negative values mean that predicted relative frequency is lower
# at higher similarity, i.e. evidence for similarity avoidance.

diff_draws <- epred_high - epred_low

# Summarize the posterior difference distribution for each
# position combination.

effect_similarity <- tibble(
  position = grid_low$position,
  Estimate = apply(diff_draws, 2, mean),
  Q2.5 = apply(diff_draws, 2, quantile, probs = 0.025),
  Q97.5 = apply(diff_draws, 2, quantile, probs = 0.975)
) %>%
  mutate(
    position_label = recode(position, !!!position_levels_pretty)
  ) %>%
  mutate(
    position_label = factor(
      position_label,
      levels = unname(position_levels_pretty)
    )
  )

p_effect <- ggplot(
  effect_similarity,
  aes(
    x = position_label,
    y = Estimate,
    ymin = Q2.5,
    ymax = Q97.5
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = 2,
    linewidth = 0.5
  ) +
  geom_pointrange() +
  labs(
    x = "Position",
    y = "Change in predicted relative frequency\n(high similarity - low similarity)",
    title = "Estimated effect of similarity by position",
    subtitle = "Negative values indicate lower predicted relative frequency at higher similarity"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title.position = "plot",
    strip.background = element_rect(fill = "grey85"),
    panel.grid.minor = element_blank()
  )

print(p_effect)

ggsave(
  filename = "arabic_offset_ratio_conditional_effect.png",
  plot = p_effect,
  width = 8,
  height = 6,
  dpi = 300
)