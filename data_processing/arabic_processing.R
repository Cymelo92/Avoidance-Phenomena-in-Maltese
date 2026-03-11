library(tidyverse)

# -----------------------------------------------------------------------------
# Load source data
# -----------------------------------------------------------------------------
wehr <- read_tsv("raw_data/wehr-arabic.csv")

# For using the same processing functions
source("data_processing/maltese_processing.R")

# -----------------------------------------------------------------------------
# Helper: convert Arabic triliteral roots to IPA
# -----------------------------------------------------------------------------
# This function:
# 1. normalizes some Arabic orthographic variants,
# 2. removes diacritics,
# 3. maps letters to IPA,
# 4. inserts hyphens between consonants for root-style output.
arabic_root_to_ipa <- function(x) {
  # Normalize orthographic variants
  x <- x %>%
    str_replace_all("[إأآٱ]", "ا") %>%  # normalize alef variants
    str_replace_all("ى", "ي") %>%      # alif maqṣūra -> ya
    str_replace_all("ؤ|ئ", "ء") %>%    # hamza carriers -> bare hamza
    str_replace_all("ـ", "") %>%       # remove tatweel
    str_replace_all("[\\u064B-\\u0652\\u0670]", "")  # remove diacritics
  
  # Arabic letter -> IPA mapping
  ipa_map <- c(
    "ء"="ʔ","ء"="ʔ", "ا"="ʔ",
    "ب"="b", "ت"="t", "ث"="θ",
    "ج"="ʤ", "ح"="ħ", "خ"="x",
    "د"="d", "ذ"="ð", "ر"="r",
    "ز"="z", "س"="s", "ش"="ʃ",
    "ص"="ṣ", "ض"="ḍ", "ط"="ṭ", "ظ"="ẓ",
    "ع"="ʕ", "غ"="ɣ",
    "ف"="f", "ق"="q", "ك"="k",
    "ل"="l", "م"="m", "ن"="n",
    "ه"="h", "و"="w", "ي"="j"
  )
  
  x %>%
    str_replace_all(ipa_map) %>%
    str_replace_all("(?<=.)(?=.)", "-")
}

# -----------------------------------------------------------------------------
# Extract canonical triliteral roots
# -----------------------------------------------------------------------------
# Keep only entries that:
# - are marked as roots,
# - consist of exactly 3 Arabic characters,
# - do not contain tāʾ marbūṭa or shadda,
# - are unique.
tri_roots <- wehr %>%
  filter(
    IsRoot == 1,
    str_length(ArabicWord) == 3,
    !str_detect(ArabicWord, "ة"),
    !str_detect(ArabicWord, "ّ")
  ) %>%
  distinct(ArabicWord) %>%
  mutate(root_ipa = arabic_root_to_ipa(ArabicWord))

# Separate root consonants and clean up the resulting data frame
arabic_roots <- 
  prepare_roots(tri_roots, root_ipa) |> 
  mutate(
    across(
      c(C1, C2, C3),
      ~ str_replace_all(.x, c(
        "ʤ" = "dʒ",
        "ṣ" = "sˤ",
        "ḍ" = "dˤ",
        "ṭ" = "tˤ",
        "ẓ" = "ðˤ"
      ))
    )
  ) |> 
  select(-ArabicWord) |> 
  rename_with(~ paste0(.x, "_curr"))

# -----------------------------------------------------------------------------
# Save output
# -----------------------------------------------------------------------------
write_csv(arabic_roots, "processed_data/arabic_roots.csv")
