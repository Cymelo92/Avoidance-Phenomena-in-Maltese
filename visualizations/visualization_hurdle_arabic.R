# ============================================================
# Arabic hurdle model: hurdle-component prediction plots
# ============================================================
# This script produces two figures from the HURDLE COMPONENT
# (hu) of the fitted Arabic hurdle negative binomial model:
#
# 1. A continuous prediction plot showing the model-implied
#    probability that a consonant pair is unattested (zero count)
#    across the observed similarity range.
#
# 2. A contrast plot showing the estimated change in the
#    predicted probability of non-occurrence between low
#    similarity (z = -1) and high similarity (z = 1).
#
# In the hurdle component:
#   - higher predicted values mean a pair is more likely to be absent
#   - positive high-minus-low contrasts mean that higher similarity
#     increases the probability of non-occurrence
#
# These plots therefore visualize similarity avoidance as lexical
# exclusion rather than as reduced frequency among attested pairs.
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

m_hurdle <- read_rds("results/arabic_hurdle.rds")

all_pairs <- m_hurdle$data %>%
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
# 3. Plot 1: continuous predictions for the hurdle component
# ------------------------------------------------------------
# We create a dense prediction grid over the observed similarity
# range, crossed with all positions.
#
# The hurdle component models the probability that freq = 0,
# i.e. that a consonant pair is unattested in the lexicon.

grid_hu <- expand.grid(
  similarity_z = seq(
    min(all_pairs$similarity_z, na.rm = TRUE),
    max(all_pairs$similarity_z, na.rm = TRUE),
    length.out = 200
  ),
  position = levels(all_pairs$position)
) %>%
  as_tibble()

# Obtain posterior summaries for the hurdle probability on the
# response scale. In brms, dpar = "hu" with scale = "response"
# returns the predicted probability of a zero count.

epred_hu <- fitted(
  m_hurdle,
  newdata = grid_hu,
  re_formula = NA,
  dpar = "hu",
  scale = "response",
  summary = TRUE
) %>%
  as.data.frame()

# Combine prediction grid with posterior summaries and create
# readable labels for plotting.

plot_hu <- bind_cols(grid_hu, epred_hu) %>%
  mutate(
    p_zero = Estimate,
    p_zero_low = Q2.5,
    p_zero_high = Q97.5,
    position_label = recode(position, !!!position_levels_pretty)
  ) %>%
  mutate(
    position_label = factor(
      position_label,
      levels = unname(position_levels_pretty)
    )
  )

# Continuous hurdle-probability plot.
# Grey ribbons show uncertainty.

p_hu <- ggplot(
  plot_hu,
  aes(
    x = similarity_z,
    y = p_zero,
    ymin = p_zero_low,
    ymax = p_zero_high
  )
) +
  geom_ribbon(alpha = 0.15, fill = "grey70", colour = NA) +
  geom_line() +
  facet_grid(. ~ position_label) +
  labs(
    x = "Phonetic similarity (z-score)",
    y = "Predicted probability of non-occurrence",
    title = "Predicted non-occurrence probability by phonetic similarity",
    subtitle = "Hurdle component of the hurdle negative binomial model"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title.position = "plot",
    strip.background = element_rect(fill = "grey85"),
    panel.grid.minor = element_blank()
  )

print(p_hu)

ggsave(
  filename = "arabic_hurdle_zero_probability.png",
  plot = p_hu,
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 4. Plot 2: contrast between low and high similarity
# ------------------------------------------------------------
# This plot summarizes the effect of similarity on the hurdle
# probability by comparing posterior predictions at:
#   low similarity  = z = -1
#   high similarity = z = 1
#
# Positive values mean that higher similarity increases the
# probability that a pair is unattested.

grid_low <- expand.grid(
  similarity_z = -1,
  position = levels(all_pairs$position)
) %>%
  as_tibble()

grid_high <- grid_low %>%
  mutate(similarity_z = 1)

# Get full posterior draws for the hurdle probability in both grids.
# We use posterior_epred() with dpar = "hu" so that we can compute
# draw-by-draw contrasts.

epred_hu_low <- posterior_epred(
  m_hurdle,
  newdata = grid_low,
  re_formula = NA,
  dpar = "hu"
)

epred_hu_high <- posterior_epred(
  m_hurdle,
  newdata = grid_high,
  re_formula = NA,
  dpar = "hu"
)

# Compute posterior draws for the change in predicted probability
# of non-occurrence from low to high similarity.
#
# Positive values indicate that higher similarity makes a pair more
# likely to be unattested.

diff_hu_draws <- epred_hu_high - epred_hu_low

# Summarize the posterior difference distribution for each
# position combination.

effect_hu <- tibble(
  position = grid_low$position,
  Estimate = apply(diff_hu_draws, 2, mean),
  Q2.5 = apply(diff_hu_draws, 2, quantile, probs = 0.025),
  Q97.5 = apply(diff_hu_draws, 2, quantile, probs = 0.975)
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

p_hu_effect <- ggplot(
  effect_hu,
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
    y = "Change in predicted non-occurrence probability\n(high similarity - low similarity)",
    title = "Estimated effect of similarity on non-occurrence probability",
    subtitle = "Positive values indicate higher non-occurrence probability at higher similarity"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title.position = "plot",
    strip.background = element_rect(fill = "grey85"),
    panel.grid.minor = element_blank()
  )

print(p_hu_effect)

ggsave(
  filename = "arabic_hurdle_zero_probability_effect.png",
  plot = p_hu_effect,
  width = 8,
  height = 6,
  dpi = 300
)