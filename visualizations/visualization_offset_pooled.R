# ============================================================
# Pooled offset model: prediction plots
# ============================================================
# This script produces two figures from the pooled cross-language
# zero-truncated negative binomial offset model:
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

fit_offset_pooled <- read_rds("results/offset_pooled_language_model.rds")

pair_counts2 <- fit_offset_pooled$data %>%
  mutate(
    position = factor(as.character(position), levels = c("12", "23", "13")),
    language = factor(as.character(language), levels = c("maltese", "arabic", "italian"))
  )

# ------------------------------------------------------------
# 2. Common relabeling helpers
# ------------------------------------------------------------

position_levels_pretty <- c(
  "12" = "C1-C2 (adjacent)",
  "23" = "C2-C3 (adjacent)",
  "13" = "C1-C3 (non-adjacent)"
)

language_levels_pretty <- c(
  "maltese" = "Maltese",
  "arabic" = "Arabic",
  "italian" = "Italian"
)

line_cols <- c(
  "Maltese" = "cornflowerblue",
  "Arabic" = "#D95F02",
  "Italian" = "#1B9E77"
)

# ------------------------------------------------------------
# 3. Plot 1: continuous predictions across similarity
# ------------------------------------------------------------
# We create a dense prediction grid over the observed similarity
# range, crossed with all positions and languages.
#
# expected_count is set to 1 for all rows. Because the fitted model
# contains offset(log(expected_count)), this makes the fitted mean
# directly interpretable as model-implied relative frequency,
# conditional on attestation.

grid_offset <- expand.grid(
  similarity_z = seq(
    min(pair_counts2$similarity_z, na.rm = TRUE),
    max(pair_counts2$similarity_z, na.rm = TRUE),
    length.out = 200
  ),
  position = levels(pair_counts2$position),
  language = levels(pair_counts2$language),
  expected_count = 1
) %>%
  as_tibble() %>%
  mutate(
    position = factor(position, levels = c("12", "23", "13")),
    language = factor(language, levels = c("maltese", "arabic", "italian"))
  )

# Obtain posterior fitted values on the response scale.
# Because expected_count = 1, these can be interpreted as
# predicted relative frequency among attested pairs.

epred_offset <- fitted(
  fit_offset_pooled,
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
    position_label = recode(as.character(position), !!!position_levels_pretty),
    language_label = recode(as.character(language), !!!language_levels_pretty)
  ) %>%
  mutate(
    position_label = factor(
      position_label,
      levels = unname(position_levels_pretty)
    ),
    language_label = factor(
      language_label,
      levels = unname(language_levels_pretty)
    )
  )

# Continuous prediction plot.
# Coloured ribbons and lines distinguish languages.

p_ratio <- ggplot(
  plot_ratio,
  aes(
    x = similarity_z,
    y = rel_freq,
    ymin = rel_freq_low,
    ymax = rel_freq_high,
    colour = language_label,
    fill = language_label,
    group = language_label
  )
) +
  geom_ribbon(alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_grid(. ~ position_label) +
  scale_colour_manual(values = line_cols) +
  scale_fill_manual(values = line_cols) +
  labs(
    x = "Phonetic similarity (z-score)",
    y = "Predicted relative frequency\namong attested pairs",
    colour = "Language",
    fill = "Language",
    title = "Predicted relative frequency by phonetic similarity",
    subtitle = "Pooled zero-truncated offset model. Expected count fixed at 1"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title.position = "plot",
    strip.background = element_rect(fill = "grey85"),
    panel.grid.minor = element_blank()
  )

print(p_ratio)

ggsave(
  filename = "pooled_offset_ratio_conditional.png",
  plot = p_ratio,
  width = 8.5,
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
  language = levels(pair_counts2$language),
  expected_count = 1
) %>%
  as_tibble() %>%
  mutate(
    position = factor(position, levels = c("12", "23", "13")),
    language = factor(language, levels = c("maltese", "arabic", "italian"))
  )

grid_high <- grid_low %>%
  mutate(similarity_z = 1)

# Get full posterior expected values for both grids.
# We use posterior_epred() rather than fitted() because we want
# draw-by-draw differences, not just posterior summaries.

epred_low <- posterior_epred(
  fit_offset_pooled,
  newdata = grid_low,
  re_formula = NA
)

epred_high <- posterior_epred(
  fit_offset_pooled,
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
# position-language combination.

effect_similarity <- tibble(
  position = grid_low$position,
  language = grid_low$language,
  Estimate = apply(diff_draws, 2, mean),
  Q2.5 = apply(diff_draws, 2, quantile, probs = 0.025),
  Q97.5 = apply(diff_draws, 2, quantile, probs = 0.975)
) %>%
  mutate(
    position_label = recode(as.character(position), !!!position_levels_pretty),
    language_label = recode(as.character(language), !!!language_levels_pretty)
  ) %>%
  mutate(
    position_label = factor(
      position_label,
      levels = unname(position_levels_pretty)
    ),
    language_label = factor(
      language_label,
      levels = unname(language_levels_pretty)
    )
  )

# Position dodge separates languages so that they do not overlap
# within each positional category.

pd <- position_dodge(width = 0.5)

p_effect <- ggplot(
  effect_similarity,
  aes(
    x = position_label,
    y = Estimate,
    ymin = Q2.5,
    ymax = Q97.5,
    colour = language_label
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = 2,
    linewidth = 0.5
  ) +
  geom_pointrange(position = pd) +
  scale_colour_manual(values = line_cols) +
  labs(
    x = "Position",
    y = "Change in predicted relative frequency\n(high similarity - low similarity)",
    colour = "Language",
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
  filename = "pooled_offset_effect.png",
  plot = p_effect,
  width = 8.5,
  height = 6,
  dpi = 300
)