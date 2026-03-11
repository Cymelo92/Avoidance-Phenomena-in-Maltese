library(tidyverse)

# -----------------------------------------------------------------------------
# 1. Read Italian IPA forms and normalize them
# -----------------------------------------------------------------------------
# Input: a text file where each line contains one or more IPA-transcribed words.
# The pipeline:
# - splits on whitespace
# - collapses geminates
# - removes stress and length marks
# - normalizes some segment symbols
# - drops forms containing unwanted symbols
# - removes vowels
# - removes final /r/
# - keeps unique consonant skeletons
#
# IMPORTANT:
# This script makes specific choices about how to normalize the Italian data:
# - geminates are collapsed
# - final r is removed (infinitive suffix -re)
# - problematic symbols are excluded


it_words <- read_lines("raw_data/italian_verbs_ipa.txt") %>%
  str_split("\\s+") %>%
  unlist(use.names = FALSE) %>%
  str_replace_all(c(
    "[ːˈˌ]"   = "",      # remove length and stress marks
    "ŋ"       = "n",
    "ts"      = "ʦ",
    "dz"      = "ʣ",
    "tʃ"      = "ʧ",
    "dʒ"      = "ʤ",
    "ɾ"       = "r"
  )) %>%
  str_replace_all(c(
    "(.)\\1+" = "\\1"    # collapse geminates
  )) %>%
  discard(~ str_detect(.x, "ə|ʒ|θ|h|c|ʝ|ɹ")) %>%
  discard(~ str_detect(.x, "\\(|\\)")) %>%
  str_remove_all("[aeiouAEIOUɪɛɔʊàèéìòóù̪-]") %>%  # remove vowels
  str_remove("r$") %>%                            # optional: remove final r
  unique()

# Keep only forms with at least 3 consonants, since triplets require 3 segments
it_words <- it_words %>%
  discard(~ str_length(.x) < 3)

# -----------------------------------------------------------------------------
# 2. Helper: tokenize a cleaned consonant string into segments
# -----------------------------------------------------------------------------
# Because your normalization already converts affricates to single symbols
# like ʧ, ʤ, ʦ, ʣ, simple character splitting is enough here.
#
# If later you keep digraph-like representations, this function would need
# to become more sophisticated.
split_segments <- function(x) {
  str_split(x, "", simplify = FALSE)[[1]]
}

# -----------------------------------------------------------------------------
# 3. Helper: extract all overlapping consonant triplets from a word
# -----------------------------------------------------------------------------
# Example:
# "strk" -> ["str", "trk"] as segment triplets
#
# Output is a tibble with:
# - word: original consonant skeleton
# - triple_id: order of the triplet within the word
# - C1, C2, C3: segments in the triplet
# - triple: a space-separated representation of the triplet
extract_triplets <- function(word) {
  segs <- split_segments(word)
  n <- length(segs)
  
  if (n < 3) {
    return(tibble())
  }
  
  map_dfr(seq_len(n - 2), function(i) {
    tibble(
      word = word,
      triple_id = i,
      C1 = segs[i],
      C2 = segs[i + 1],
      C3 = segs[i + 2],
      triple = paste(segs[i], segs[i + 1], segs[i + 2])
    )
  })
}

# -----------------------------------------------------------------------------
# 4. Build the Italian triplet table
# -----------------------------------------------------------------------------
# Each word may yield multiple overlapping triplets.
italian_triplets <- map_dfr(it_words, extract_triplets)

italian_triplets <- italian_triplets %>%
  mutate(
    across(
      where(is.character),
      ~ str_replace_all(.x, c(
        "ʦ" = "ts",
        "ʣ" = "dz",
        "ʧ" = "tʃ",
        "ʤ" = "dʒ"
      ))
    )
  ) |> 
  rename_with(~ paste0(.x, "_curr"), .cols = c(C1, C2, C3)) |> 
  select(C1_curr, C2_curr, C3_curr)

write_csv(italian_triplets, "processed_data/italian_roots.csv")

