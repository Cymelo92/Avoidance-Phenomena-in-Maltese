# ============================================================
# Visualizing the Maltese offset model in brms
# Model object name assumed here: fit_maltese_offset
#
# This script gives several complementary visualizations for the
# offset negative-binomial model:
#
#   count ~ similarity_z * position * origin + offset(log(expected_count))
#
# The key interpretive point is:
# this model is NOT a raw-count model. Because of the offset,
# it asks whether a consonant pair occurs more or less often than
# expected under the independence baseline.
#
# So the most useful plots are not just "counts by similarity".
# The most useful plots are:
#
#   1. Predicted count for a standardized expected_count = 1
#      -> this is effectively a relative count plot
#
#   2. Model-implied observed/expected ratio
#      -> values below 1 mean underattestation relative to expectation
#      -> values above 1 mean overattestation relative to expectation
#
#   3. Raw observed/expected ratios in the data
#      -> descriptive, not model-based
#
#   4. Raw observed/expected ratios with the posterior model fit overlaid
#      -> lets you compare the noisy data to the model trend
#
#   5. A compact "effect of similarity" plot by position and origin
#      -> shows how much predicted counts change when similarity increases
#
#   6. Standard conditional_effects() plots
#      -> quick diagnostic / presentation plots from brms
#
#   7. Posterior predictive checks
#      -> checks whether the model reproduces the count distribution well
#
# Assumptions:
# - your fitted model is called: fit_maltese_offset
# - your modeling data frame is called: pair_counts2
# - pair_counts2 contains:
#     count, expected_count, similarity_z, position, origin
# ============================================================

library(brms)
library(tidyverse)
library(posterior)

# ------------------------------------------------------------
# 0. Sanity checks
# ------------------------------------------------------------
# Make sure the main variables are in the expected form.
# position and origin should be factors with stable level order.
# This matters for faceting and for correct reference categories.

fit_maltese_offset <- read_rds("results/maltese_offset.rds")

pair_counts2 <- fit_maltese_offset$data

pair_counts2 <- pair_counts2 %>%
  mutate(
    position = factor(position, levels = c("12", "23", "13")),
    origin = factor(origin, levels = c("sem", "loan"))
  )

# ------------------------------------------------------------
# 1. Build a prediction grid
# ------------------------------------------------------------
# We create a grid of similarity values across the observed range,
# crossed with all positions and origins.
#
# IMPORTANT:
# We set expected_count = 1 in the prediction grid.
# Why?
# Because the offset is multiplicative. Setting expected_count = 1
# makes the predicted count directly interpretable as a relative rate:
# "how often would the pair occur, relative to expectation?"
#
# If expected_count = 1, then:
#   predicted count = predicted observed/expected ratio

grid_offset <- expand.grid(
  similarity_z = seq(
    min(pair_counts2$similarity_z, na.rm = TRUE),
    max(pair_counts2$similarity_z, na.rm = TRUE),
    length.out = 100
  ),
  position = levels(pair_counts2$position),
  origin = levels(pair_counts2$origin),
  expected_count = 1
) %>%
  as_tibble()

# ------------------------------------------------------------
# 2. Posterior expected counts on the response scale
# ------------------------------------------------------------
# fitted(..., scale = "response") returns posterior summaries for the
# expected observed count under the model.
#
# Because expected_count = 1 in the grid, these are easy to interpret
# as relative counts / observed-to-expected tendencies.

epred_offset <- fitted(
  fit_maltese_offset,
  newdata = grid_offset,
  re_formula = NA,    # population-level effects only; exclude random effects if there were any
  scale = "response",
  summary = TRUE
) %>%
  as.data.frame()

plot_pred <- bind_cols(grid_offset, epred_offset)

# ------------------------------------------------------------
# 3. Plot 1: Predicted count by similarity
# ------------------------------------------------------------
# Interpretation:
# This shows the model's predicted count for a hypothetical pair whose
# expected_count is fixed at 1.
#
# Because expected_count = 1, this plot is already close to an
# observed/expected interpretation:
#   - values below 1 suggest underattestation relative to expectation
#   - values above 1 suggest overattestation relative to expectation
#
# This is a clean "main effects + interactions" visualization.

p_pred_count <- ggplot(
  plot_pred,
  aes(x = similarity_z, y = Estimate, ymin = Q2.5, ymax = Q97.5)
) +
  geom_ribbon(alpha = 0.20) +
  geom_line(linewidth = 0.9) +
  facet_grid(origin ~ position) +
  labs(
    x = "Similarity (z-score)",
    y = "Predicted count (expected_count = 1)",
    title = "Predicted pair count by similarity",
    subtitle = "Offset negative-binomial model"
  ) +
  theme_bw()

print(p_pred_count)

# ------------------------------------------------------------
# 4. Plot 2: Model-implied observed/expected ratio
# ------------------------------------------------------------
# This is often the single most interpretable plot for the offset model.
#
# Since the model contains an offset(log(expected_count)), the quantity
# people usually care about is:
#
#   predicted observed count / expected count
#
# Since expected_count = 1 in the prediction grid, this is numerically
# identical to the predicted count. But plotting it explicitly as a ratio
# makes the interpretation much clearer.
#
# The dashed horizontal line at 1 is crucial:
#   - ratio < 1  -> underattested relative to expectation
#   - ratio > 1  -> overattested relative to expectation

plot_ratio <- plot_pred %>%
  mutate(
    oe_ratio = Estimate / expected_count,
    oe_ratio_low = Q2.5 / expected_count,
    oe_ratio_high = Q97.5 / expected_count
  )

p_ratio <- ggplot(
  plot_ratio,
  aes(x = similarity_z, y = oe_ratio, ymin = oe_ratio_low, ymax = oe_ratio_high)
) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_ribbon(alpha = 0.20) +
  geom_line(linewidth = 0.9) +
  facet_grid(origin ~ position) +
  labs(
    x = "Similarity (z-score)",
    y = "Observed / expected ratio",
    title = "Model-implied observed/expected ratio by similarity",
    subtitle = "Values below 1 indicate underattestation"
  ) +
  theme_bw()

print(p_ratio)

# ------------------------------------------------------------
# 5. Prepare raw observed/expected ratios in the data
# ------------------------------------------------------------
# This is the descriptive, non-model-based version.
#
# It is useful because it lets you see how noisy the actual data are.
# But it should not replace the posterior model plots.
#
# Raw O/E values can be unstable, especially when expected counts are small.

pair_counts2 <- pair_counts2 %>%
  mutate(
    oe_ratio = count / expected_count
  )

# ------------------------------------------------------------
# 6. Plot 3: Raw observed/expected ratios by similarity
# ------------------------------------------------------------
# Interpretation:
# This is a descriptive plot of the actual data.
# It shows where individual pair-position-origin cells fall.
#
# This plot is useful to show the cloud of observations.
# It is not a formal inferential summary.
#
# The dashed line at 1 again marks exact match to expectation.

p_raw_oe <- ggplot(
  pair_counts2,
  aes(x = similarity_z, y = oe_ratio)
) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_point(alpha = 0.40) +
  facet_grid(origin ~ position) +
  labs(
    x = "Similarity (z-score)",
    y = "Observed / expected ratio",
    title = "Raw observed/expected ratios by similarity",
    subtitle = "Descriptive plot of the data"
  ) +
  theme_bw()

print(p_raw_oe)

# ------------------------------------------------------------
# 7. Plot 4: Raw O/E ratios with posterior model fit overlaid
# ------------------------------------------------------------
# This is often the best combined figure:
# it shows both the noisy data and the posterior trend.
#
# Interpretation:
# - points = actual pair cells
# - ribbon/line = model-estimated trend with uncertainty
# - horizontal line at 1 = exact expectation
#
# This lets you show that the model is smoothing over noisy raw ratios.

p_raw_plus_model <- ggplot(
  pair_counts2,
  aes(x = similarity_z, y = oe_ratio)
) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_point(alpha = 0.35) +
  geom_ribbon(
    data = plot_ratio,
    aes(
      x = similarity_z,
      y = oe_ratio,
      ymin = oe_ratio_low,
      ymax = oe_ratio_high,
      group = interaction(origin, position)
    ),
    inherit.aes = FALSE,
    alpha = 0.20
  ) +
  geom_line(
    data = plot_ratio,
    aes(
      x = similarity_z,
      y = oe_ratio,
      group = interaction(origin, position)
    ),
    inherit.aes = FALSE,
    linewidth = 0.9
  ) +
  facet_grid(origin ~ position) +
  labs(
    x = "Similarity (z-score)",
    y = "Observed / expected ratio",
    title = "Observed/expected ratios with posterior model fit",
    subtitle = "Points are raw data; ribbon and line are model predictions"
  ) +
  theme_bw()

print(p_raw_plus_model)

# ------------------------------------------------------------
# 8. Plot 5: Compact similarity-effect plot by position and origin
# ------------------------------------------------------------
# Because the model contains interactions, the coefficient table is hard
# to read directly. A clearer summary is to compare predictions at
# low similarity and high similarity.
#
# Here we compare:
#   similarity_z = -1  versus  similarity_z = +1
#
# Since similarity is z-scored, this is a 2-SD contrast.
#
# Interpretation:
#   negative values -> higher similarity lowers predicted count
#   positive values -> higher similarity raises predicted count
#
# Since expected_count = 1, this is again interpretable on a relative scale.

grid_low <- expand.grid(
  similarity_z = -1,
  position = levels(pair_counts2$position),
  origin = levels(pair_counts2$origin),
  expected_count = 1
) %>%
  as_tibble()

grid_high <- grid_low %>%
  mutate(similarity_z = 1)

epred_low <- posterior_epred(
  fit_maltese_offset,
  newdata = grid_low,
  re_formula = NA
)

epred_high <- posterior_epred(
  fit_maltese_offset,
  newdata = grid_high,
  re_formula = NA
)

diff_draws <- epred_high - epred_low

effect_similarity <- tibble(
  position = grid_low$position,
  origin = grid_low$origin,
  Estimate = apply(diff_draws, 2, mean),
  Q2.5 = apply(diff_draws, 2, quantile, probs = 0.025),
  Q97.5 = apply(diff_draws, 2, quantile, probs = 0.975)
)

p_effect <- ggplot(
  effect_similarity,
  aes(x = position, y = Estimate, ymin = Q2.5, ymax = Q97.5)
) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_pointrange() +
  facet_wrap(~ origin) +
  labs(
    x = "Position",
    y = "Change in predicted count\n(high similarity - low similarity)",
    title = "Estimated effect of similarity by position and origin",
    subtitle = "Negative values indicate fewer counts at higher similarity"
  ) +
  theme_bw()

print(p_effect)

# ------------------------------------------------------------
# 9. Plot 6: Standard brms conditional_effects() plot
# ------------------------------------------------------------
# This is the quick, built-in visualization.
# It is useful for exploratory work and draft plots.
#
# Here we condition on origin and display the similarity × position
# interaction. Because this is a 3-way model, conditions = ... is needed.
#
# method = "posterior_epred" ensures we are plotting expected counts on the
# response scale rather than the linear predictor.

cond_offset <- data.frame(origin = c("sem", "loan"))

ce_offset <- conditional_effects(
  fit_maltese_offset,
  effects = "similarity_z:position",
  conditions = cond_offset,
  method = "posterior_epred",
  re_formula = NA
)

plot(ce_offset)

# ------------------------------------------------------------
# 10. Optional alternative conditional_effects() layout
# ------------------------------------------------------------
# This version conditions on position and plots similarity × origin.
# Sometimes this layout is easier to read if the origin contrast is
# the main substantive question.

cond_offset_pos <- data.frame(position = c("12", "23", "13"))

ce_offset_alt <- conditional_effects(
  fit_maltese_offset,
  effects = "similarity_z:origin",
  conditions = cond_offset_pos,
  method = "posterior_epred",
  re_formula = NA
)

plot(ce_offset_alt)

# ------------------------------------------------------------
# 11. Posterior predictive checks
# ------------------------------------------------------------
# These do not directly visualize the similarity effect.
# Instead, they check whether the model reproduces the overall count
# distribution reasonably well.
#
# Recommended:
# - default pp_check(): quick overview
# - bars: useful for discrete count data
# - rootogram: often very good for count models

pp_check(fit_maltese_offset)
pp_check(fit_maltese_offset, type = "bars", ndraws = 100)
pp_check(fit_maltese_offset, type = "rootogram", ndraws = 100)

# ------------------------------------------------------------
# 12. Optional: save the main plots
# ------------------------------------------------------------
# Uncomment if you want files written to disk.

# ggsave("maltese_offset_predicted_count.png", p_pred_count, width = 8, height = 5, dpi = 300)
# ggsave("maltese_offset_oe_ratio.png", p_ratio, width = 8, height = 5, dpi = 300)
# ggsave("maltese_offset_raw_oe.png", p_raw_oe, width = 8, height = 5, dpi = 300)
# ggsave("maltese_offset_raw_plus_model.png", p_raw_plus_model, width = 8, height = 5, dpi = 300)
# ggsave("maltese_offset_similarity_effect.png", p_effect, width = 8, height = 5, dpi = 300)

# ============================================================
# Which plots are most worth keeping?
#
# If you only want 2 figures:
#   1. p_ratio
#      -> clearest visualization of underattestation relative to expectation
#   2. p_effect
#      -> clearest summary of where the similarity effect is strongest
#
# If you want 3 figures:
#   add p_raw_plus_model
#      -> best combination of raw data + posterior fit
#
# Interpretation reminders:
#
# p_pred_count / p_ratio:
#   shows how predicted relative frequency changes with similarity
#
# p_raw_oe:
#   descriptive only; raw cells can be noisy
#
# p_raw_plus_model:
#   lets you compare noisy data with posterior trend
#
# p_effect:
#   summarizes the effect of moving from low to high similarity
#
# pp_check():
#   checks whether the model captures the count distribution at all
# ============================================================

