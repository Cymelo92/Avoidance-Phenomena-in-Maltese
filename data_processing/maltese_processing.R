library(tidyverse)

# Read source data
df <- read_csv("raw_data/verb-mt.csv")

# -----------------------------------------------------------------------------
# Helper: split a triliteral root into three consonant slots
# -----------------------------------------------------------------------------
# Expects a column containing roots formatted like "k-t-b".
# Returns the input data frame with columns C1, C2, and C3.
prepare_roots <- function(data, root_col) {
  root_col <- rlang::ensym(root_col)
  
  data %>%
    separate(!!root_col, into = c("C1", "C2", "C3"), sep = "-", remove = TRUE) %>%
    filter(
      C1 != "", C2 != "", C3 != "",
      C2 != C3
    ) %>%
    mutate(
      across(c(C1, C2, C3), ~ str_replace_all(.x, "’", "għ"))
    )
}

# -----------------------------------------------------------------------------
# Helper: map orthographic consonants to historical and current pronunciations
# -----------------------------------------------------------------------------
# Input must contain C1, C2, and C3.
# Output includes, for each consonant position:
#   - orthographic form
#   - historical form
#   - current form
map_consonants <- function(
    roots,
    orth = c("b", "d", "w", "ċ", "f", "ġ", "g", "għ", "ħ", "h", "j", "k",
             "l", "m", "n", "p", "q", "r", "s", "t", "v", "x", "ż", "z"),
    hist = c("b", "d", "w", "tʃ", "f", "dʒ", "g", "ʕ", "ħ", "h", "j", "k",
             "l", "m", "n", "p", "q", "r", "s", "t", "v", "ʃ", "z", "ts"),
    curr = c("b", "d", "w", "tʃ", "f", "dʒ", "g", "ː", "h", "0", "j", "k",
             "l", "m", "n", "p", "ʔ", "r", "s", "t", "v", "ʃ", "z", "ts")
) {
  needed <- c("C1", "C2", "C3")
  
  if (!all(needed %in% names(roots))) {
    stop("Input data frame must contain columns: C1, C2, C3")
  }
  
  if (!(length(orth) == length(hist) && length(orth) == length(curr))) {
    stop("orth, hist, and curr must have the same length")
  }
  
  map_hist <- setNames(hist, orth)
  map_curr <- setNames(curr, orth)
  
  roots %>%
    mutate(
      across(
        all_of(needed),
        list(
          hist = ~ recode(.x, !!!map_hist, .default = .x),
          curr = ~ recode(.x, !!!map_curr, .default = .x)
        ),
        .names = "{.col}_{.fn}"
      )
    ) %>%
    rename(
      C1_orth = C1,
      C2_orth = C2,
      C3_orth = C3
    ) %>%
    select(
      C1_orth, C1_hist, C1_curr,
      C2_orth, C2_hist, C2_curr,
      C3_orth, C3_hist, C3_curr
    )
}

# -----------------------------------------------------------------------------
# Extract and process Semitic triliteral roots
# -----------------------------------------------------------------------------
# Keep only distinct triliteral roots and exclude entries with separators or
# special symbols that indicate non-canonical forms.
sem_roots <- df %>%
  filter(għerq_kategorija == "trilitteru") %>%
  select(għerq) %>%
  distinct() %>%
  filter(!str_detect(għerq, "/|•|\\\\")) %>%
  prepare_roots(għerq) %>%
  map_consonants()

# -----------------------------------------------------------------------------
# Extract and process loanword-derived roots
# -----------------------------------------------------------------------------
# For rows without an annotated root:
#   1. simplify geminates in elicitation_form,
#   2. remove vowels,
#   3. insert hyphens between remaining consonants,
#   4. keep only forms with exactly three consonants.
loan_roots <- df %>%
  filter(is.na(għerq)) %>%
  mutate(
    elicitation_form = str_replace_all(elicitation_form, "(.)\\1", "\\1"),
    għerq = elicitation_form %>%
      str_remove_all(regex("[aeiou]", ignore_case = TRUE)) %>%
      str_split("") %>%
      map_chr(~ paste(.x[.x != ""], collapse = "-"))
  ) %>%
  filter(str_count(għerq, "-") == 2) %>%
  select(għerq) %>%
  prepare_roots(għerq) %>%
  map_consonants()

# -----------------------------------------------------------------------------
# Save processed data
# -----------------------------------------------------------------------------

write_csv(sem_roots, "processed_data/maltese_sem_roots.csv")
write_csv(loan_roots, "processed_data/maltese_loan_roots.csv")
