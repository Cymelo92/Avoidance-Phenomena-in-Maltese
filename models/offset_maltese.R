library(tidyverse)
library(brms)

# -----------------------------------------------------------------------------
# 1. Load and combine the Maltese datasets
# -----------------------------------------------------------------------------
# Read the two datasets used for the offset model:
# - loanword-derived roots
# - Semitic roots
#
# Add a variable marking the origin of each observation, then combine them into
# a single data frame.
ml <- read_csv2("processed_data/maltese-loan_offset-data.csv") %>%
  mutate(origin = "loan")

ms <- read_csv2("processed_data/maltese-sem_offset-data.csv") %>%
  mutate(origin = "sem")

df <- bind_rows(ml, ms)

# -----------------------------------------------------------------------------
# 2. Basic cleaning and type conversion
# -----------------------------------------------------------------------------
# Keep only the variables needed for modeling and convert them to the correct
# types:
# - root and consonants as character strings
# - position as a factor with a fixed level order
# - similarity as numeric
# - origin as a factor with "sem" as the reference level
#
# Remove any rows with missing values in these core fields.
dat_raw <- df %>%
  transmute(
    root = as.character(root),
    consonants = as.character(consonants),
    position = factor(as.character(position), levels = c("12", "23", "13")),
    similarity = as.numeric(similarity),
    origin = factor(as.character(origin), levels = c("sem", "loan"))
  ) %>%
  filter(
    !is.na(root),
    !is.na(consonants),
    !is.na(position),
    !is.na(similarity),
    !is.na(origin)
  )

# -----------------------------------------------------------------------------
# 3. Helper functions for parsing roots and consonant pairs
# -----------------------------------------------------------------------------
# The data assume roots and consonant pairs are written with separators,
# typically "-" or spaces.
#
# Examples:
# - root: "k-t-b" or "k t b"
# - pair: "k-t" or "k t"

# Split a consonant pair into its two segments.
split_pair <- function(x) {
  x <- str_trim(x)
  
  if (str_detect(x, "-")) {
    parts <- str_split(x, "-", simplify = TRUE)
  } else if (str_detect(x, "\\s+")) {
    parts <- str_split(x, "\\s+", simplify = TRUE)
  } else {
    stop(
      "Could not parse consonant pair: ", x,
      ". Expected a separator like '-' or space."
    )
  }
  
  parts <- parts[parts != ""]
  
  if (length(parts) != 2) {
    stop("Pair does not have exactly 2 consonants: ", x)
  }
  
  parts
}

# Split a triliteral root into C1, C2, and C3.
split_root <- function(x) {
  x <- str_trim(x)
  
  if (str_detect(x, "-")) {
    parts <- str_split(x, "-", simplify = TRUE)
  } else if (str_detect(x, "\\s+")) {
    parts <- str_split(x, "\\s+", simplify = TRUE)
  } else {
    stop(
      "Could not parse root: ", x,
      ". Expected root segments to be separated by '-' or space."
    )
  }
  
  parts <- parts[parts != ""]
  
  if (length(parts) != 3) {
    stop("Root does not have exactly 3 consonants: ", x)
  }
  
  tibble(C1 = parts[1], C2 = parts[2], C3 = parts[3])
}

# Given a position label, return the corresponding root slots.
# For example:
# - "12" means the pair should match C1 and C2
# - "23" means C2 and C3
# - "13" means C1 and C3
pair_slots_from_position <- function(position) {
  if (position == "12") return(c("C1", "C2"))
  if (position == "23") return(c("C2", "C3"))
  if (position == "13") return(c("C1", "C3"))
  stop("Unknown position: ", position)
}

# -----------------------------------------------------------------------------
# 4. Parse roots into C1, C2, and C3
# -----------------------------------------------------------------------------
# Parse each distinct root once, then join the parsed consonant slots back onto
# the full dataset.
root_parsed <- dat_raw %>%
  distinct(root, origin) %>%
  rowwise() %>%
  mutate(root_parts = list(split_root(root))) %>%
  unnest(root_parts) %>%
  ungroup()

dat <- dat_raw %>%
  left_join(root_parsed, by = c("root", "origin"))

# -----------------------------------------------------------------------------
# 5. Parse observed consonant pairs and verify consistency with the root
# -----------------------------------------------------------------------------
# For each row:
# - split the observed pair into its two consonants
# - determine which root slots that pair is supposed to represent, based on
#   the position code
# - check that the observed pair matches the relevant consonants in the root
pair_checked <- dat %>%
  rowwise() %>%
  mutate(pair_parts = list(split_pair(consonants))) %>%
  mutate(
    pair_c1 = pair_parts[1],
    pair_c2 = pair_parts[2]
  ) %>%
  ungroup() %>%
  select(-pair_parts) %>%
  mutate(
    root_slot1 = case_when(
      position == "12" ~ C1,
      position == "23" ~ C2,
      position == "13" ~ C1
    ),
    root_slot2 = case_when(
      position == "12" ~ C2,
      position == "23" ~ C3,
      position == "13" ~ C3
    ),
    matches_root = (pair_c1 == root_slot1 & pair_c2 == root_slot2)
  )

# Stop if any pair does not match the root-position coding.
if (any(!pair_checked$matches_root)) {
  bad_rows <- pair_checked %>%
    filter(!matches_root) %>%
    select(root, consonants, position, pair_c1, pair_c2, C1, C2, C3) %>%
    head(10)
  
  print(bad_rows)
  stop("Some consonant pairs do not match the root and position coding.")
}

dat <- pair_checked

# -----------------------------------------------------------------------------
# 6. Collapse the data to counts by pair, position, and origin
# -----------------------------------------------------------------------------
# The model is fit to aggregated counts rather than individual rows.
# For each consonant pair within each position and origin:
# - count how many times it occurs
# - carry forward its similarity value
pair_counts <- dat %>%
  group_by(consonants, position, origin) %>%
  summarise(
    count = n(),
    similarity = first(similarity),
    .groups = "drop"
  )

# Verify that similarity is unique within each aggregated cell.
sim_check <- dat %>%
  group_by(consonants, position, origin) %>%
  summarise(n_similarity = n_distinct(similarity), .groups = "drop") %>%
  filter(n_similarity > 1)

if (nrow(sim_check) > 0) {
  print(sim_check)
  stop("Similarity is not unique within some consonants × position × origin cells.")
}

# -----------------------------------------------------------------------------
# 7. Compute expected counts under an independence baseline
# -----------------------------------------------------------------------------
# The goal is to estimate how often each pair would be expected to occur if the
# two consonants combined independently, given slot-specific segment frequencies.
#
# Example:
# For position 12 in origin = sem,
# expected_count = N_pairs_(12, sem) * P(C1 = seg1 | sem) * P(C2 = seg2 | sem)

# Convert parsed roots to long format so slot-specific segment frequencies can
# be estimated.
slot_freq <- root_parsed %>%
  pivot_longer(
    cols = c(C1, C2, C3),
    names_to = "slot",
    values_to = "segment"
  ) %>%
  count(origin, slot, segment, name = "slot_segment_n") %>%
  group_by(origin, slot) %>%
  mutate(slot_segment_p = slot_segment_n / sum(slot_segment_n)) %>%
  ungroup()

# Count the total number of observed pair tokens in each origin × position stratum.
position_totals <- dat %>%
  count(origin, position, name = "N_position_origin")

# Split each consonant pair into its two members and determine which slots they
# correspond to.
pair_counts2 <- pair_counts %>%
  rowwise() %>%
  mutate(pair_parts = list(split_pair(consonants))) %>%
  mutate(
    seg1 = pair_parts[1],
    seg2 = pair_parts[2]
  ) %>%
  ungroup() %>%
  select(-pair_parts) %>%
  mutate(
    slot1 = case_when(
      position == "12" ~ "C1",
      position == "23" ~ "C2",
      position == "13" ~ "C1"
    ),
    slot2 = case_when(
      position == "12" ~ "C2",
      position == "23" ~ "C3",
      position == "13" ~ "C3"
    )
  ) %>%
  left_join(
    slot_freq %>%
      rename(slot1 = slot, seg1 = segment, p1 = slot_segment_p),
    by = c("origin", "slot1", "seg1")
  ) %>%
  left_join(
    slot_freq %>%
      rename(slot2 = slot, seg2 = segment, p2 = slot_segment_p),
    by = c("origin", "slot2", "seg2")
  ) %>%
  left_join(position_totals, by = c("origin", "position")) %>%
  mutate(expected_count = N_position_origin * p1 * p2)

# Safety checks for the expected counts:
# - they must not be missing
# - they must be positive, because the model uses log(expected_count) as an offset
if (any(is.na(pair_counts2$expected_count))) {
  bad <- pair_counts2 %>%
    filter(is.na(expected_count))
  print(bad)
  stop("Some expected counts are NA. Check root parsing and segment coding.")
}

if (any(pair_counts2$expected_count <= 0)) {
  stop("Some expected counts are <= 0, which makes log(offset) invalid.")
}

# -----------------------------------------------------------------------------
# 8. Standardize similarity
# -----------------------------------------------------------------------------
# Standardizing similarity makes the coefficient easier to interpret and usually
# helps model fitting.
pair_counts2 <- pair_counts2 %>%
  mutate(similarity_z = as.numeric(scale(similarity)))

# -----------------------------------------------------------------------------
# 9. Define priors
# -----------------------------------------------------------------------------
# The model predicts pair counts using:
# - similarity_z
# - position
# - origin
# - all interactions among them
#
# The offset adjusts the expected count under independence, so the model asks
# whether pairs occur more or less often than that baseline predicts.
priors <- c(
  prior(normal(0, 1), class = "b"),
  prior(normal(0, 1), class = "Intercept"),
  prior(exponential(1), class = "shape")
)

# -----------------------------------------------------------------------------
# 10. Fit Bayesian negative binomial model with offset
# -----------------------------------------------------------------------------
# Response:
# - count: observed frequency of each consonant pair
#
# Predictor structure:
# - similarity_z * position * origin
#
# Offset:
# - log(expected_count), which encodes the independence baseline
#
# Truncation:
# - lb = 1 because the aggregated table includes only observed pairs
fit <- brm(
  formula = bf(
    count | trunc(lb = 1) ~ similarity_z * position * origin + offset(log(expected_count))
  ),
  data = pair_counts2,
  family = negbinomial(),
  prior = priors,
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  seed = 1234
)

# -----------------------------------------------------------------------------
# 11. Save fitted model
# -----------------------------------------------------------------------------
write_rds(fit, "results/maltese_offset.rds")