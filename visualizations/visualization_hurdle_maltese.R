# ============================================================
# Visualizing the Maltese hurdle model in brms
# Model object name assumed here: m_hurdle
#
# This script gives several complementary visualizations for the
# hurdle negative-binomial model:
#
#   freq ~ similarity_z * position * origin + (1 | C1) + (1 | C2)
#   hu   ~ similarity_z * position * origin + (1 | C1) + (1 | C2)
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
#      pair-position-origin combination is unattested.
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
#   4. Effect-of-similarity plots by position and origin
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
#     similarity_z, position, origin, freq, C1, C2
# ============================================================

library(brms)
library(tidyverse)
library(posterior)

# ------------------------------------------------------------
# 0. Sanity checks
# ------------------------------------------------------------
# Make sure the core grouping variables are factors with the intended
# level order. This matters for faceting and for interpretation of
# reference categories in the model.

m_hurdle <- read_rds("results/maltese_hurdle.rds")

all_pairs <- m_hurdle$data

all_pairs <- all_pairs %>%
  mutate(
    position = factor(position, levels = c("12", "23", "13")),
    origin = factor(origin, levels = c("sem", "loan"))
  )

# ------------------------------------------------------------
# 1. Build a prediction grid
# ------------------------------------------------------------
# We create a grid across the observed range of similarity_z,
# crossed with all positions and origins.
#
# These are the combinations for which we will ask the model:
# "What does it predict here?"

grid <- expand.grid(
  similarity_z = seq(
    min(all_pairs$similarity_z, na.rm = TRUE),
    max(all_pairs$similarity_z, na.rm = TRUE),
    length.out = 100
  ),
  position = levels(all_pairs$position),
  origin = levels(all_pairs$origin)
) %>%
  as_tibble()

# ------------------------------------------------------------
# 2. Plot 1: Overall expected frequency
# ------------------------------------------------------------
# This is the most general plot for the hurdle model.
#
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
#
# This is the broadest summary because it collapses both mechanisms
# into a single expected count.

epred <- fitted(
  m_hurdle,
  newdata = grid,
  re_formula = NA,   # population-level effect only; omit random effects
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
  facet_grid(origin ~ position) +
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
# This is one of the most important plots for your analysis.
#
# The hurdle parameter hu is the probability that the count is zero.
# In your application, that means:
#
#   probability the pair-position-origin combination is unattested
#
# So this plot answers:
#   "As similarity increases, how does the probability that a pair is
#    absent from the lexicon change?"
#
# Interpretation:
#   higher values = more likely to be unattested
#   lower values  = more likely to occur
#
# If similarity increases this curve, that is direct evidence that
# similarity pushes combinations toward lexical absence.

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
  facet_grid(origin ~ position) +
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
# The previous plot (p_mean) combines both parts of the hurdle model.
# But sometimes we want to isolate the count mechanism alone.
#
# In a hurdle model, the count part models the expected frequency
# conditional on the count being positive.
#
# That is:
#   E[y | y > 0]
#
# In brms, this corresponds to the count-model mean parameter mu.
# So plotting mu shows what happens to frequency among attested pairs.
#
# This plot answers:
#   "Given that a pair does occur, how does similarity affect how
#    frequent it is?"
#
# This is distinct from the hurdle plot:
#   - p_hu is about whether the pair occurs at all
#   - p_mu is about how frequent it is once attested

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
  facet_grid(origin ~ position) +
  labs(
    x = "Similarity (z-score)",
    y = "Mean positive frequency",
    title = "Positive-count mean by similarity",
    subtitle = "Expected frequency conditional on attestation"
  ) +
  theme_bw()

print(p_mu)

# ------------------------------------------------------------
# 5. Plot 4: Group-specific effect of similarity on expected frequency
# ------------------------------------------------------------
# This plot is a compact summary of the interaction structure.
#
# Because your model includes similarity * position * origin, the
# coefficient table is hard to read directly. A more interpretable
# approach is to compare the model's predictions at:
#
#   similarity_z = -1   and   similarity_z = +1
#
# Since similarity_z is standardized, this is a 2-standard-deviation
# contrast from low to high similarity.
#
# We compute the posterior difference:
#
#   expected frequency at high similarity
#   minus
#   expected frequency at low similarity
#
# Interpretation:
#   negative values -> higher similarity lowers expected frequency
#   positive values -> higher similarity raises expected frequency
#
# This plot summarizes the net effect of similarity for each
# position × origin cell.

grid_lo <- expand.grid(
  similarity_z = -1,
  position = levels(all_pairs$position),
  origin = levels(all_pairs$origin)
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
  origin = grid_lo$origin,
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
  facet_wrap(~ origin) +
  labs(
    x = "Position",
    y = "Change in expected frequency\n(high similarity - low similarity)",
    title = "Effect of similarity on expected frequency",
    subtitle = "Negative values indicate lower expected frequency at higher similarity"
  ) +
  theme_bw()

print(p_effect_mean)

# ------------------------------------------------------------
# 6. Plot 5: Group-specific effect of similarity on unattestedness
# ------------------------------------------------------------
# This is often the most theoretically important compact plot.
#
# We do the same low-vs-high similarity comparison as above,
# but now for the hurdle probability hu.
#
# So the posterior difference is:
#
#   probability unattested at high similarity
#   minus
#   probability unattested at low similarity
#
# Interpretation:
#   positive values -> higher similarity increases unattestedness
#   negative values -> higher similarity decreases unattestedness
#
# This plot is especially useful because your coefficient table
# suggested that similarity affects non-occurrence more strongly
# than positive frequencies.

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
  origin = grid_lo$origin,
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
  facet_wrap(~ origin) +
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
# These are built-in brms visualizations. They are useful for quick
# inspection and draft figures.
#
# Because your model includes a 3-way interaction, you need
# conditions = ... to fix one variable and display the other two.
#
# Here, we condition on origin and plot the similarity × position
# interaction.

cond <- data.frame(origin = c("sem", "loan"))

# Overall expected count on the response scale
ce_mean <- conditional_effects(
  m_hurdle,
  effects = "similarity_z:position",
  conditions = cond,
  method = "posterior_epred",
  re_formula = NA
)

plot(ce_mean)

# Hurdle probability, i.e. probability of being unattested
ce_hu <- conditional_effects(
  m_hurdle,
  effects = "similarity_z:position",
  conditions = cond,
  dpar = "hu",
  method = "posterior_epred",
  re_formula = NA
)

plot(ce_hu)

# ------------------------------------------------------------
# 8. Optional alternative conditional_effects() layout
# ------------------------------------------------------------
# Instead of conditioning on origin, you can condition on position and
# look at similarity × origin.
#
# This layout is often useful if the origin contrast is central to the
# theoretical question.

cond2 <- data.frame(position = c("12", "23", "13"))

ce_mean2 <- conditional_effects(
  m_hurdle,
  effects = "similarity_z:origin",
  conditions = cond2,
  method = "posterior_epred",
  re_formula = NA
)

plot(ce_mean2)

ce_hu2 <- conditional_effects(
  m_hurdle,
  effects = "similarity_z:origin",
  conditions = cond2,
  dpar = "hu",
  method = "posterior_epred",
  re_formula = NA
)

plot(ce_hu2)

# ------------------------------------------------------------
# 9. Posterior predictive checks
# ------------------------------------------------------------
# These plots do not show the similarity effect directly.
# Instead, they ask whether the model reproduces the overall shape of
# the observed data well.
#
# For the hurdle model, these are especially important because they
# tell you whether the model is getting:
#   - the large number of zeros
#   - the low counts
#   - the right tail of the count distribution
#
# Recommended:
#   - default pp_check(): general overview
#   - bars: useful for discrete counts
#   - rootogram: very useful for count distributions with many zeros

pp_check(m_hurdle)
pp_check(m_hurdle, type = "bars", ndraws = 100)
pp_check(m_hurdle, type = "rootogram", ndraws = 100)

# ------------------------------------------------------------
# 10. Optional: save the main plots
# ------------------------------------------------------------
# Uncomment if you want image files written to disk.

# ggsave("maltese_hurdle_expected_frequency.png", p_mean, width = 8, height = 5, dpi = 300)
# ggsave("maltese_hurdle_probability_unattested.png", p_hu, width = 8, height = 5, dpi = 300)
# ggsave("maltese_hurdle_positive_count_mean.png", p_mu, width = 8, height = 5, dpi = 300)
# ggsave("maltese_hurdle_effect_similarity_frequency.png", p_effect_mean, width = 8, height = 5, dpi = 300)
# ggsave("maltese_hurdle_effect_similarity_unattestedness.png", p_effect_hu, width = 8, height = 5, dpi = 300)

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
#   compact summary of how similarity changes expected frequency by group
#
# p_effect_hu:
#   compact summary of how similarity changes unattestedness by group
#
# pp_check():
#   checks whether the hurdle model captures zeros and count distribution
# ============================================================