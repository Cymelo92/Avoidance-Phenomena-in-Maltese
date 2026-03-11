# ============================================================
# Visualizing the Arabic hurdle model in brms
# Model object name assumed here: m_hurdle
#
# This script gives several complementary visualizations for the
# hurdle negative-binomial model:
#
#   freq ~ similarity_z * position + (1 | C1) + (1 | C2)
#   hu   ~ similarity_z * position + (1 | C1) + (1 | C2)
#
# IMPORTANT CONCEPTUAL POINT:
# A hurdle model has TWO linked parts:
#
#   1. Count part (mu):
#      This models the size of the count, conditional on the count
#      being positive.
#
#   2. Hurdle part (hu):
#      This models whether the observation is zero at all.
#      In this application, that corresponds to whether a consonant
#      pair-position combination is unattested.
#
# So, unlike a simple count model, the hurdle model lets you separate:
#
#   - Does similarity make a pair more likely to be absent?
#   - If the pair does occur, does similarity make it less frequent?
#
# The most useful plots are therefore:
#
#   1. Overall expected frequency on the response scale
#      -> combines both the hurdle and count parts
#
#   2. Hurdle probability
#      -> probability that a pair is unattested (i.e. frequency = 0)
#
#   3. Positive-count mean
#      -> expected count conditional on the pair being attested
#
#   4. Effect-of-similarity plots by position
#      -> compact summaries of where the similarity effect is strongest
#
#   5. Standard conditional_effects() plots
#      -> quick built-in brms visualizations
#
#   6. Posterior predictive checks
#      -> checks whether the model reproduces the zero inflation and
#         count distribution reasonably well
#
# Assumptions:
# - your fitted hurdle model is called: m_hurdle
# - your modeling data frame is called: all_pairs
# - all_pairs contains:
#     similarity_z, position, freq, C1, C2
# ============================================================

library(brms)
library(tidyverse)
library(posterior)

# ------------------------------------------------------------
# 0. Sanity checks
# ------------------------------------------------------------
# Make sure position is a factor with the intended level order.
# This matters for plotting and for interpretation of reference
# categories in the model.

m_hurdle <- read_rds("results/arabic_hurdle.rds")

all_pairs <- m_hurdle$data

all_pairs <- all_pairs %>%
  mutate(
    position = factor(position, levels = c("12", "23", "13"))
  )

# ------------------------------------------------------------
# 1. Build a prediction grid
# ------------------------------------------------------------
# Create a grid across the observed range of similarity_z,
# crossed with all positions.
#
# These are the combinations for which we will ask the model:
# "What does it predict here?"

grid <- expand.grid(
  similarity_z = seq(
    min(all_pairs$similarity_z, na.rm = TRUE),
    max(all_pairs$similarity_z, na.rm = TRUE),
    length.out = 100
  ),
  position = levels(all_pairs$position)
) %>%
  as_tibble()

# ------------------------------------------------------------
# 2. Plot 1: Overall expected frequency
# ------------------------------------------------------------
# fitted(..., scale = "response") gives the posterior expected value
# of the response on the observed-data scale.
#
# For a hurdle model, this combines:
#   - probability of being non-zero
#   - expected count when non-zero
#
# So this plot answers:
#   "As similarity increases, what happens to the expected observed
#    frequency overall?"

epred <- fitted(
  m_hurdle,
  newdata = grid,
  re_formula = NA,
  scale = "response",
  summary = TRUE
) %>%
  as.data.frame()

plot_mean <- bind_cols(grid, epred)

p_mean <- ggplot(
  plot_mean,
  aes(x = similarity_z, y = Estimate, ymin = Q2.5, ymax = Q97.5)
) +
  geom_ribbon(alpha = 0.20) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ position) +
  labs(
    x = "Similarity (z-score)",
    y = "Expected frequency",
    title = "Overall expected frequency by similarity",
    subtitle = "Hurdle negative-binomial model; population-level effects only"
  ) +
  theme_bw()

print(p_mean)

# ------------------------------------------------------------
# 3. Plot 2: Hurdle probability = probability unattested
# ------------------------------------------------------------
# The hurdle parameter hu is the probability that the count is zero.
# In this application, that means:
#
#   probability the pair-position combination is unattested
#
# So this plot answers:
#   "As similarity increases, how does the probability that a pair is
#    absent from the lexicon change?"
#
# Interpretation:
#   higher values = more likely to be unattested
#   lower values  = more likely to occur

hu_pred <- fitted(
  m_hurdle,
  newdata = grid,
  re_formula = NA,
  dpar = "hu",
  scale = "response",
  summary = TRUE
) %>%
  as.data.frame()

plot_hu <- bind_cols(grid, hu_pred)

p_hu <- ggplot(
  plot_hu,
  aes(x = similarity_z, y = Estimate, ymin = Q2.5, ymax = Q97.5)
) +
  geom_ribbon(alpha = 0.20) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ position) +
  labs(
    x = "Similarity (z-score)",
    y = "Probability unattested",
    title = "Hurdle probability by similarity",
    subtitle = "Higher values mean the pair is more likely to be absent"
  ) +
  theme_bw()

print(p_hu)

# ------------------------------------------------------------
# 4. Plot 3: Positive-count mean
# ------------------------------------------------------------
# This isolates the count mechanism alone.
#
# In a hurdle model, the count part models the expected frequency
# conditional on the count being positive:
#
#   E[y | y > 0]
#
# In brms, this corresponds to the count-model mean parameter mu.
#
# This plot answers:
#   "Given that a pair does occur, how does similarity affect how
#    frequent it is?"

mu_draws <- posterior_linpred(
  m_hurdle,
  newdata = grid,
  re_formula = NA,
  dpar = "mu",
  transform = TRUE
)

mu_summary <- tibble(
  Estimate = apply(mu_draws, 2, mean),
  Q2.5 = apply(mu_draws, 2, quantile, probs = 0.025),
  Q97.5 = apply(mu_draws, 2, quantile, probs = 0.975)
)

plot_mu <- bind_cols(grid, mu_summary)

p_mu <- ggplot(
  plot_mu,
  aes(x = similarity_z, y = Estimate, ymin = Q2.5, ymax = Q97.5)
) +
  geom_ribbon(alpha = 0.20) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ position) +
  labs(
    x = "Similarity (z-score)",
    y = "Mean positive frequency",
    title = "Positive-count mean by similarity",
    subtitle = "Expected frequency conditional on attestation"
  ) +
  theme_bw()

print(p_mu)

# ------------------------------------------------------------
# 5. Plot 4: Position-specific effect of similarity on expected frequency
# ------------------------------------------------------------
# Because the model includes similarity * position, the coefficient
# table is harder to interpret directly than a contrast plot.
#
# We compare the model's predictions at:
#   similarity_z = -1   and   similarity_z = +1
#
# Since similarity_z is standardized, this is a 2-SD contrast from
# low to high similarity.
#
# Interpretation:
#   negative values -> higher similarity lowers expected frequency
#   positive values -> higher similarity raises expected frequency

grid_lo <- expand.grid(
  similarity_z = -1,
  position = levels(all_pairs$position)
) %>%
  as_tibble()

grid_hi <- grid_lo %>%
  mutate(similarity_z = 1)

ep_lo <- posterior_epred(
  m_hurdle,
  newdata = grid_lo,
  re_formula = NA
)

ep_hi <- posterior_epred(
  m_hurdle,
  newdata = grid_hi,
  re_formula = NA
)

diff_ep <- ep_hi - ep_lo

effect_mean <- tibble(
  position = grid_lo$position,
  Estimate = apply(diff_ep, 2, mean),
  Q2.5 = apply(diff_ep, 2, quantile, probs = 0.025),
  Q97.5 = apply(diff_ep, 2, quantile, probs = 0.975)
)

p_effect_mean <- ggplot(
  effect_mean,
  aes(x = position, y = Estimate, ymin = Q2.5, ymax = Q97.5)
) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_pointrange() +
  labs(
    x = "Position",
    y = "Change in expected frequency\n(high similarity - low similarity)",
    title = "Effect of similarity on expected frequency",
    subtitle = "Negative values indicate lower expected frequency at higher similarity"
  ) +
  theme_bw()

print(p_effect_mean)

# ------------------------------------------------------------
# 6. Plot 5: Position-specific effect of similarity on unattestedness
# ------------------------------------------------------------
# This is the same low-vs-high similarity comparison as above,
# but now for the hurdle probability hu.
#
# The posterior difference is:
#
#   probability unattested at high similarity
#   minus
#   probability unattested at low similarity
#
# Interpretation:
#   positive values -> higher similarity increases unattestedness
#   negative values -> higher similarity decreases unattestedness

hu_lo <- posterior_epred(
  m_hurdle,
  newdata = grid_lo,
  re_formula = NA,
  dpar = "hu"
)

hu_hi <- posterior_epred(
  m_hurdle,
  newdata = grid_hi,
  re_formula = NA,
  dpar = "hu"
)

diff_hu <- hu_hi - hu_lo

effect_hu <- tibble(
  position = grid_lo$position,
  Estimate = apply(diff_hu, 2, mean),
  Q2.5 = apply(diff_hu, 2, quantile, probs = 0.025),
  Q97.5 = apply(diff_hu, 2, quantile, probs = 0.975)
)

p_effect_hu <- ggplot(
  effect_hu,
  aes(x = position, y = Estimate, ymin = Q2.5, ymax = Q97.5)
) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_pointrange() +
  labs(
    x = "Position",
    y = "Change in hurdle probability\n(high similarity - low similarity)",
    title = "Effect of similarity on unattestedness",
    subtitle = "Positive values indicate higher unattestedness at higher similarity"
  ) +
  theme_bw()

print(p_effect_hu)

# ------------------------------------------------------------
# 7. Standard brms conditional_effects() plots
# ------------------------------------------------------------
# These are built-in brms visualizations for quick inspection.
#
# Here we plot the similarity × position interaction.

ce_mean <- conditional_effects(
  m_hurdle,
  effects = "similarity_z:position",
  method = "posterior_epred",
  re_formula = NA
)

plot(ce_mean)

# Hurdle probability, i.e. probability of being unattested
ce_hu <- conditional_effects(
  m_hurdle,
  effects = "similarity_z:position",
  dpar = "hu",
  method = "posterior_epred",
  re_formula = NA
)

plot(ce_hu)

# ------------------------------------------------------------
# 8. Posterior predictive checks
# ------------------------------------------------------------
# These do not show the similarity effect directly.
# They evaluate whether the model reproduces the overall shape of
# the observed data, especially:
#   - the number of zeros
#   - the low counts
#   - the right tail of the count distribution

pp_check(m_hurdle)
pp_check(m_hurdle, type = "bars", ndraws = 100)
pp_check(m_hurdle, type = "rootogram", ndraws = 100)

# ------------------------------------------------------------
# 9. Optional: save the main plots
# ------------------------------------------------------------
# Uncomment if you want image files written to disk.

# ggsave("arabic_hurdle_expected_frequency.png", p_mean, width = 8, height = 5, dpi = 300)
# ggsave("arabic_hurdle_probability_unattested.png", p_hu, width = 8, height = 5, dpi = 300)
# ggsave("arabic_hurdle_positive_count_mean.png", p_mu, width = 8, height = 5, dpi = 300)
# ggsave("arabic_hurdle_effect_similarity_frequency.png", p_effect_mean, width = 8, height = 5, dpi = 300)
# ggsave("arabic_hurdle_effect_similarity_unattestedness.png", p_effect_hu, width = 8, height = 5, dpi = 300)

# ============================================================
# Which plots are most worth keeping?
#
# If you only want 2 figures:
#   1. p_mean
#      -> broadest summary: what similarity does to overall expected frequency
#   2. p_hu
#      -> most theoretically important: what similarity does to unattestedness
#
# If you want 3 figures:
#   add p_effect_hu
#      -> compact summary of where similarity most strongly increases absence
#
# Interpretation reminders:
#
# p_mean:
#   overall expected observed frequency (combines occurrence and count size)
#
# p_hu:
#   probability that a pair is absent from the lexicon
#
# p_mu:
#   frequency among attested pairs only
#
# p_effect_mean:
#   compact summary of how similarity changes expected frequency by position
#
# p_effect_hu:
#   compact summary of how similarity changes unattestedness by position
#
# pp_check():
#   checks whether the hurdle model captures zeros and count distribution
# ============================================================