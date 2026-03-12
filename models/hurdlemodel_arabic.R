library(tidyverse)
library(brms)

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
# Read the dataset used for the Arabic hurdle model.
all_pairs <- read_csv2("processed_data/arabic_hurdle-data.csv")

# -----------------------------------------------------------------------------
# Prepare variables for modeling
# -----------------------------------------------------------------------------
# - position is treated as an ordered factor with a fixed reference order.
# - C1 and C2 are grouping factors for random intercepts.
# - similarity is kept numeric and also standardized for modeling.
# - freq and occurs are cast to integer to match count/binary data structure.
all_pairs <- all_pairs %>%
  mutate(
    position = factor(position, levels = c(12, 23, 13)),
    C1 = factor(C1),
    C2 = factor(C2),
    similarity = as.numeric(similarity),
    similarity_z = as.numeric(scale(similarity)),
    freq = as.integer(freq),
    occurs = as.integer(occurs)
  )

# -----------------------------------------------------------------------------
# Fit Bayesian hurdle negative binomial model
# -----------------------------------------------------------------------------
# The model has two parts:
# 1. Count model: predicts non-zero frequency values.
# 2. Hurdle model (hu): predicts whether a zero vs non-zero outcome occurs.
#
# Fixed effects:
# - similarity_z
# - position
# - their interaction
#
# Random effects:
# - varying intercepts for C1 and C2 in both model components

# Set priors
priors <- c(
  prior(normal(0, 1), class = "b"),
  prior(normal(0, 1.5), class = "Intercept"),
  prior(normal(0, 1), class = "b", dpar = "hu"),
  prior(normal(0, 1.5), class = "Intercept", dpar = "hu"),
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "shape")
)

m_hurdle <- brm(
  bf(
    freq ~ similarity_z * position + (1 | C1) + (1 | C2),
    hu ~ similarity_z * position + (1 | C1) + (1 | C2)
  ),
  data = all_pairs,
  family = hurdle_negbinomial(),
  prior = priors,
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  seed = 123
)

# -----------------------------------------------------------------------------
# Save fitted model
# -----------------------------------------------------------------------------
# Store the model object for later inspection, prediction, or plotting.
write_rds(m_hurdle, "results/arabic_hurdle.rds")