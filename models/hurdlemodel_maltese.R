library(tidyverse)
library(brms)

# -----------------------------------------------------------------------------
# Load and combine data
# -----------------------------------------------------------------------------
# Read the two Maltese datasets and label each row by origin:
# - "loan" for loanword-derived roots
# - "sem" for Semitic roots
ml <- read_csv2("processed_data/maltese-loan_hurdle-data.csv") %>%
  mutate(origin = "loan")

ms <- read_csv2("processed_data/maltese-sem_hurdle-data.csv") %>%
  mutate(origin = "sem")

all_pairs <- bind_rows(ml, ms)

# -----------------------------------------------------------------------------
# Prepare variables for modeling
# -----------------------------------------------------------------------------
# - position is treated as a factor with a fixed level order.
# - origin distinguishes Semitic vs loan roots.
# - C1 and C2 are grouping factors for random intercepts.
# - similarity is kept numeric and also standardized.
# - freq and occurs are cast to integers for count/binary modeling.
all_pairs <- all_pairs %>%
  mutate(
    position = factor(position, levels = c(12, 23, 13)),
    origin = factor(origin, levels = c("sem", "loan")),
    C1 = factor(C1),
    C2 = factor(C2),
    similarity = as.numeric(similarity),
    similarity_z = as.numeric(scale(similarity)),
    freq = as.integer(freq),
    occurs = as.integer(occurs)
  )

# -----------------------------------------------------------------------------
# Set priors
# -----------------------------------------------------------------------------
# Priors are specified separately for readability.
# The same general prior structure is used for:
# - the count component
# - the hurdle (hu) component
priors <- c(
  prior(normal(0, 1), class = "b"),
  prior(normal(0, 1.5), class = "Intercept"),
  prior(normal(0, 1), class = "b", dpar = "hu"),
  prior(normal(0, 1.5), class = "Intercept", dpar = "hu"),
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "shape")
)

# -----------------------------------------------------------------------------
# Fit Bayesian hurdle negative binomial model
# -----------------------------------------------------------------------------
# This model differs from the Arabic model in including origin as an additional
# predictor, along with all interactions among:
# - similarity_z
# - position
# - origin
#
# As before, both the count and hurdle components include random intercepts for
# C1 and C2.
m_hurdle <- brm(
  bf(
    freq ~ similarity_z * position * origin + (1 | C1) + (1 | C2),
    hu   ~ similarity_z * position * origin + (1 | C1) + (1 | C2)
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
# Save the model object for later summaries, diagnostics, and plotting.
write_rds(m_hurdle, "results/maltese_hurdle.rds")