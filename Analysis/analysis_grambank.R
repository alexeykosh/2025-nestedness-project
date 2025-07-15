#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Script for Grambank analysis ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This code allows for the visualization of Grambank data using several plots 
# to show the distribution of features. It also merges the data with Glottolog 
# in order to retrieve the names of the language families.

# REQUIRED INPUT FILES:
# 1. Grambank data (statistical curated dataset in CLDF format from the article):
#    - values.csv => Grammatical feature values per language
#    - languages.csv => Metadata for Grambank languages (incl. level-1 family)
#
# 2. Glottolog data:
#    - languoid.csv => Contains mappings from family codes to full family names

## ==== Libraries ====
library(ggplot2)
library(dplyr)

## ==== Load Grambank data ====
df_values <- read.csv("values.csv", sep = ",", header = TRUE, 
                      stringsAsFactors = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Clean and re-factor the column 'value' ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ==== Replace non-0/1/2/?/NaN values with "other" ====
valids <- c("0", "1", "2", "?", "NaN")
df_values$value_clean <- with(df_values,
                              ifelse(is.na(value) | value %in% valids,
                                     as.character(value),
                                     "other")
)

## ==== Convert to factor ====
df_values$value_clean <- factor(
  df_values$value_clean,
  levels = c("0", "1", "2", "?", "other"),
  exclude = NULL
)

## ==== Barplot of cleaned column 'value' ====
counts <- table(df_values$value_clean, useNA = "ifany")
labels <- names(counts)
labels[is.na(labels) | labels == ""] <- "NA"

barplot(
  counts,
  names.arg = labels,
  main      = "Barplot of the column 'value'",
  xlab      = "Value category",
  ylab      = "Count",
  las       = 2,
  col       = "cadetblue",
  border    = "white",
  cex.names = 0.8
)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Proportion of NA and '?' per feature ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ==== Compute stats by feature ====
stats_feat <- df_values %>%
  group_by(new.name) %>%
  summarise(
    total    = n(),
    na_count = sum(is.na(value)),
    qm_count = sum(value == "?"),
    prop_na  = na_count / total,
    prop_qm  = qm_count / total
  ) %>%
  ungroup()

## ==== Histogram of NA proportions ====
hist(
  stats_feat$prop_na,
  breaks = seq(0, 1, by = 0.05),
  main   = "Proportion of NA by Feature",
  xlab   = "Proportion of NA",
  ylab   = "Number of Features",
  col    = "cadetblue",
  border = "white",
  xlim   = c(0, 1)
)

## ==== Histogram of '?' proportions ====
hist(
  stats_feat$prop_qm,
  breaks = seq(0, 1, by = 0.05),
  main   = "Proportion of '?' by Feature",
  xlab   = "Proportion of '?'",
  ylab   = "Number of Features",
  col    = "cadetblue",
  border = "white",
  xlim   = c(0, 1)
)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Identify features with only {0,1,?,NA} ("valid" features) ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ==== Features with only valid values ====
valids_simple <- c("0", "1", "?")
feature_vals <- df_values %>%
  mutate(value = as.character(value)) %>%
  group_by(new.name) %>%
  summarise(vals = list(unique(na.omit(value)))) %>%
  ungroup() %>%
  mutate(is_simple = sapply(vals, function(v) all(v %in% valids_simple)))

n_simple <- sum(feature_vals$is_simple)
cat("Number of features with only 0, 1, ?, or NA:", n_simple, "\n")

invalid_features <- feature_vals$new.name[!feature_vals$is_simple]
cat("Features not in that group:\n")
print(invalid_features)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Distribution of languages informed per feature ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ==== Count distinct languages per feature (non-NA, non-?) ====
lang_per_feature <- df_values %>%
  filter(!is.na(value), value != "?") %>%
  group_by(new.name) %>%
  summarise(n_lang = n_distinct(glottocode)) %>%
  ungroup()

## ==== Histogram of languages per feature ====
hist(
  lang_per_feature$n_lang,
  breaks = seq(0, max(lang_per_feature$n_lang, na.rm = TRUE) + 50, by = 50),
  main  = "Distribution of number of languages per feature (without  ? or NA)",
  xlab  = "Number of languages with valid value (not ? and NA)",
  ylab = "Number of features",
  col  = "cadetblue",
  border = "white"
)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Distribution of fully-informed features per language ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ==== Count features with no NA or '?' per language ====
lang_complete <- df_values %>%
  filter(!is.na(value), value != "?") %>%
  group_by(glottocode) %>%
  summarise(n_features = n_distinct(new.name)) %>%
  ungroup()

cat("Sample of complete features per language:\n")
print(head(lang_complete))

## ==== Histogram of complete features per language ====
hist(
  lang_complete$n_features,
  breaks = "Sturges",
  main = "Number of fully-informed (no ? or NA) features per language",
  xlab  = "Count of features without any NA or '?'",
  ylab = "Number of languages",
  col   = "cadetblue",
  border = "white"
)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Number of languages per level-1 family (full names) ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ==== 2. Load the Grambank languages.csv to get level1 family codes ====
languages <- read.csv("languages.csv",
                      sep = ",",
                      header = TRUE,
                      stringsAsFactors = FALSE)

## ==== 3. Keep only glottocodes present in values.csv ====
langs_in_values <- unique(df_values$glottocode)
langs_families <- languages %>%
  filter(glottocode %in% langs_in_values) %>%
  select(glottocode, family_code = level1)

## ==== 4. Count languages per family_code ====
family_counts <- langs_families %>%
  distinct(glottocode, family_code) %>%  
  group_by(family_code) %>%
  summarise(n_lang = n()) %>%
  filter(n_lang >= 4) %>%   # keep only families with at least 4 languages
  arrange(n_lang)                   

print(family_counts)


## ==== Load the Glottolog languoid.csv to map family_code = full name ====
languoid <- read.csv("languoid.csv",
                     sep = ",",      # c’est un CSV classique
                     header = TRUE,
                     stringsAsFactors = FALSE)

## ==== Join to retrieve full family name ====
family_named <- family_counts %>%
  left_join(
    languoid %>% select(id, full_name = name),
    by = c("family_code" = "id")
  )

print(family_named)

## ==== Plot barplot with full_name ====
ggplot(family_named, aes(x = reorder(full_name, n_lang), y = n_lang)) +
  geom_col(fill = "cadetblue", color = "white") +
  coord_flip() +
  labs(
    title = "Number of languages by level-1 family (>= 4 languages)",
    x = "Language family",
    y  = "Number of languages"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(hjust = 1),
    plot.title  = element_text(hjust = 0.5)
  )
