# Libraries
library(dplyr)       
library(ggplot2)    
library(tidyr)       
library(vegan)       

theme_set(theme_bw())

# Setting main parameters
ALPHA_ <- 0.005
N_ITER_ <- 10

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# I. ONE GLOBAL MATRIX ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## 1.1 Load data ----
df_values <- read.csv(
  "data/cldf-datasets-phoible-f36deac/cldf/values.csv",
  sep = ",", header = TRUE
)

## 1.2 Select a single inventory per language ----
# Apply the same logic as in nested_test, but globally
selected_combinations <- df_values %>%
  group_by(Language_ID, Contribution_ID) %>%
  summarise(
    phoneme_list = list(unique(Value)),
    n_phonemes   = n_distinct(Value),
    .groups      = "drop"
  ) %>%
  group_by(Language_ID) %>%
  group_modify(~ {
    if (nrow(.x) == 1) {
      .x
    } else if (nrow(.x) == 2) {
      # for two inventories: choose the one with fewer phonemes
      min_count <- min(.x$n_phonemes)
      candidates <- filter(.x, n_phonemes == min_count)
      slice_sample(candidates, n = 1)
    } else {
      # for more than two inventories: maximize overlap, then minimize phoneme count
      overlaps <- sapply(seq_len(nrow(.x)), function(i) {
        sum(sapply(seq_len(nrow(.x)), function(j) {
          if (i != j) length(intersect(
            .x$phoneme_list[[i]], .x$phoneme_list[[j]]
          )) else 0
        }))
      })
      cand1 <- .x[overlaps == max(overlaps), ]
      min_count <- min(cand1$n_phonemes)
      cand2 <- filter(cand1, n_phonemes == min_count)
      slice_sample(cand2, n = 1)
    }
  }) %>%
  ungroup() %>%
  select(Language_ID, Contribution_ID)

## 1.3 Build the global matrix ----
df_selected <- df_values %>%
  inner_join(
    selected_combinations,
    by = c("Language_ID", "Contribution_ID")
  )

df_one_hot <- df_selected %>%
  select(Language_ID, Value) %>%
  mutate(count = 1) %>%
  pivot_wider(
    names_from  = Value,
    values_from = count,
    values_fill = list(count = 0)
  )

mat_global <- as.matrix(df_one_hot[, -1])
rownames(mat_global) <- df_one_hot$Language_ID

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# II. RUN NESTEDNESS TESTS ON GLOBAL MATRIX ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## 2.1 NODF tests ----
res_nodf_r00 <- oecosimu(
  mat_global,
  nestfun     = nestednodf,
  method      = "r00",
  nsimul      = N_ITER_,
  parallel    = -1,
  alternative = "two.sided"
)

res_nodf_c0 <- oecosimu(
  mat_global,
  nestfun     = nestednodf,
  method      = "c0",
  nsimul      = N_ITER_,
  parallel    = -1,
  alternative = "two.sided"
)

## 2.2 Temperature tests ----
res_temp_r00 <- oecosimu(
  mat_global,
  nestfun     = nestedtemp,
  method      = "r00",
  nsimul      = N_ITER_,
  parallel    = -1,
  alternative = "two.sided"
)

res_temp_c0 <- oecosimu(
  mat_global,
  nestfun     = nestedtemp,
  method      = "c0",
  nsimul      = N_ITER_,
  parallel    = -1,
  alternative = "two.sided"
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# III. RESULTS AS DATAFRAMES ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Helper to build a results dataframe from an oecosimu object
build_df <- function(res, measure, stat_index) {
  sims   <- res$oecosimu$simulated
  pval   <- res$oecosimu$pval[stat_index]
  real_v <- unname(res$statistic$statistic[stat_index])
  data.frame(
    Measure = measure,
    Type    = c(rep("simulated", length(sims)), "real"),
    Value   = c(sims, real_v),
    p_value = c(rep(pval, length(sims)), pval)
  )
}

## 3.1 NODF results ----
df_nodf_r00 <- build_df(res_nodf_r00, "NODF", 2)
df_nodf_c0  <- build_df(res_nodf_c0,  "NODF", 2)

## 3.2 Temperature results ----
df_temp_r00 <- build_df(res_temp_r00, "Temperature", 1)
df_temp_c0  <- build_df(res_temp_c0,  "Temperature", 1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IV. GAUSSIAN DISTRIBUTION PLOT ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Function to plot the global Gaussian distributions
plot_distribution_global <- function(df_r00, df_c0,
                                     measure_label, x_label, x_lim) {
  sim_r00 <- df_r00 %>% filter(Type == "simulated") %>% pull(Value)
  sim_c0  <- df_c0  %>% filter(Type == "simulated") %>% pull(Value)
  mean_r00 <- mean(sim_r00); sd_r00 <- sd(sim_r00)
  mean_c0  <- mean(sim_c0);  sd_c0  <- sd(sim_c0)
  real_r00 <- df_r00 %>% filter(Type == "real") %>% pull(Value)
  real_c0  <- df_c0  %>% filter(Type == "real") %>% pull(Value)
  
  x_seq <- seq(x_lim[1], x_lim[2], length.out = 200)
  df_density <- bind_rows(
    data.frame(x = x_seq,
               y = dnorm(x_seq, mean_r00, sd_r00),
               baseline = "r00"),
    data.frame(x = x_seq,
               y = dnorm(x_seq, mean_c0, sd_c0),
               baseline = "c0")
  )
  
  ggplot(df_density, aes(x, y, fill = baseline, color = baseline)) +
    geom_area(alpha = 0.3, position = "identity") +
    geom_line(size = 1) +
    # add real‐value vertical lines
    geom_vline(
      xintercept = real_r00,
      linetype   = "dashed",
      size       = 1.1,
      color      = "red"
    ) +
    geom_vline(
      xintercept = real_c0,
      linetype   = "dashed",
      size       = 1.1,
      color      = "red"
    ) +
    scale_fill_manual(
      values = c("r00" = "#56B", "c0" = "orange"),
      guide  = FALSE
    ) +
    scale_color_manual(
      name   = "Baseline",
      values = c("r00" = "#56B", "c0" = "orange")
    ) +
    labs(
      x     = x_label,
      y     = "Density"
    ) +
    scale_x_continuous(limits = x_lim) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
}

## 4.1 Plot and save NODF distribution ----
dist_nodf_global <- plot_distribution_global(
  df_nodf_r00, df_nodf_c0,
  measure_label = "NODF",
  x_label       = "NODF",
  x_lim         = range(df_nodf_r00$Value[df_nodf_r00$Type == "simulated"])
)
ggsave("dist_nodf_global.png", dist_nodf_global,
       width = 8, height = 6, bg = "white")

## 4.2 Plot and save Temperature distribution ----
dist_temp_global <- plot_distribution_global(
  df_temp_r00, df_temp_c0,
  measure_label = "Temperature",
  x_label       = "Temperature",
  x_lim         = range(df_temp_r00$Value[df_temp_r00$Type == "simulated"])
)
ggsave("dist_temp_global.png", dist_temp_global,
       width = 8, height = 6, bg = "white")
