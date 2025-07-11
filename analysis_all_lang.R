#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ==== Script for the analysis of the nestedness of phonological ====
# inventories across the whole phoible dataset as one single matrix 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This code performs the same analysis as the nestedness analysis of phonological 
# inventories by language families, but this time using the entire PHOIBLE 
# dataset combined into a single matrix.

# ==== Libraries ====
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)

# ==== Setting main parameters ====
# Number of iterations
N_ITER_ <- 1000
# set parallel options to the computer's number of cores minus 2
options(mc.cores = max(1, parallel::detectCores() - 2))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# I. NESTEDNESS METRICS ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# In this section, we create the global matrix with all phoible languages.
# We then calculate the degree of nestedness for this matrix using NODF 
# and Temperature from the oecosimu package. 
# Finally, we use baselines and simulations (r00 and c0) to compare our results.

# This section can be run in one go and will compute both Temperature and NODF

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Useful functions ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Function to select one inventory per language (same logic as the familiy analysis)
select_inventories <- function(df_values) {
  ## Step 1: Select one inventory per language 
  selected_combinations <- df_values %>%
    group_by(Language_ID, Contribution_ID) %>%
    summarise(
      phoneme_list = list(unique(Value)),
      n_phonemes = n_distinct(Value),
      .groups = 'drop'
    ) %>%
    group_by(Language_ID) %>%
    group_modify(~ {
      if (nrow(.x) == 1) {
        ## Only one inventory: keep it
        return(.x)
      } else if (nrow(.x) == 2) {
        ## For 2 inventories: select the one with fewer phonemes 
        ## random tie-break if equal
        min_count <- min(.x$n_phonemes)
        candidates <- .x %>% filter(n_phonemes == min_count)
        return(candidates %>% slice_sample(n = 1))
      } else {
        ## For more than 2 inventories: calculate the degree of overlap
        ## If equal, select the inventory with fewer phonemes
        ## random tie-break if still equal
        overlaps <- sapply(1:nrow(.x), function(i) {
          sum(sapply(1:nrow(.x), function(j) {
            if (i != j) {
              length(intersect(.x$phoneme_list[[i]], .x$phoneme_list[[j]]))
            } else {
              0
            }
          }))
        })
        max_overlap <- max(overlaps)
        candidates <- .x[overlaps == max_overlap, ]
        min_count <- min(.x$n_phonemes)
        candidates <- .x %>% filter(n_phonemes == min_count)
        return(candidates %>% slice_sample(n = 1))
      }
    }) %>%
    ungroup() %>%
    select(Language_ID, Contribution_ID)
  ## Step 2: Retrieve matching rows
  df_selected <- df_values %>%
    dplyr::inner_join(selected_combinations,
                      by = c("Language_ID", "Contribution_ID"))
  ## Step 3: One-hot encode phonemes
  df_one_hot <- df_selected %>%
    dplyr::select(Language_ID, Value) %>%
    dplyr::mutate(count = 1) %>%
    tidyr::pivot_wider(
      names_from   = Value,
      values_from  = count,
      values_fill  = list(count = 0)
    )
  ## Convert to matrix format
  df_one_hot_matrix <- as.matrix(df_one_hot[, -1])
  rownames(df_one_hot_matrix) <- df_one_hot$Language_ID
  ## Compute dimensions + fill
  nrow_ <- nrow(df_one_hot_matrix)
  ncol_ <- ncol(df_one_hot_matrix)
  fill  <- sum(df_one_hot_matrix) / (nrow_ * ncol_)
  ## Return everything we need
  list(
    selected_combinations = selected_combinations,
    mat = df_one_hot_matrix,
    nrow = nrow_,
    ncol = ncol_,
    fill = fill
  )
}

## Function for nestedness (NODF) testing on the full dataset
nested_test_nodf_full <- function(df_values, shuffling_type = 'r00') {
  ## Step 1: select inventories AND build matrix
  prep   <- select_inventories(df_values)
  mat    <- prep$mat
  n_langs <- prep$nrow
  ## Step 2: simulate
  res <- oecosimu(mat, nestfun = nestednodf, method = shuffling_type,
                  nsimul = N_ITER_, parallel = TRUE)
  stats <- res$statistic$statistic   # 3 values
  pvals <- res$oecosimu$pval[3]    
  ## Step 3: assemble real + simulated
  # real
  df_real <- data.frame(
    Measure = 'NODF',
    n_langs = n_langs,
    Baseline = shuffling_type,
    Type = "real",
    NODF_columns_Value = stats[1],
    NODF_rows_Value= stats[2],
    NODF_global_Value = stats[3],
    NODF_global_p_value = pvals,
    stringsAsFactors = TRUE
  )
  # simulated
  sim_mat <- res$oecosimu$simulated
  df_sim <- as.data.frame(t(sim_mat))
  colnames(df_sim) <- c("NODF_columns_Value","NODF_rows_Value","NODF_global_Value")
  df_sim <- df_sim %>%
    mutate(iter = row_number()) %>%
    transmute(
      Measure = 'NODF',
      n_langs = n_langs,
      Baseline = shuffling_type,
      Type = "simulated",
      NODF_columns_Value,
      NODF_rows_Value,
      NODF_global_Value,
      NODF_global_p_value = pvals
    )
  # print(res)
  # print(sim_mat)
  return(bind_rows(df_real, df_sim))
}

## Function for Temperature testing on the full dataset
nested_test_temp_full <- function(df_values, shuffling_type = 'r00') {
  ## Step 1: select inventories AND build matrix
  prep   <- select_inventories(df_values)
  mat    <- prep$mat
  n_langs <- prep$nrow
  ## Step 2: simulate
  res <- oecosimu(mat, nestfun = nestedtemp, method = shuffling_type,
                  nsimul = N_ITER_, parallel = TRUE)
  stat_t <- res$statistic$statistic[1]
  p_t  <- res$oecosimu$pval[1]
  ## Step 3: assemble real + simulated
  sim <- res$oecosimu$simulated
  # real
  df_real <- data.frame(
    Measure = 'Temperature',
    n_langs = nrow(mat),
    Baseline = shuffling_type,
    Type = rep('simulated', length(sim)),
    Value = as.numeric(sim),
    p_value = as.numeric(res$oecosimu$pval[1])
  )
  # simulated
  df_sim <- data.frame(
    Measure = 'Temperature',
    n_langs = nrow(mat),
    Baseline = shuffling_type,
    Type = 'real',
    Value = as.numeric(res$statistic$statistic[1]),
    p_value = as.numeric(res$oecosimu$pval[1])
  )
  # print(res)
  # print(sim)
  return(bind_rows(df_real, df_sim))
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Loading and processing data ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Read Phoible data
df_values <- read.csv("values.csv",
                      sep = ",", header = TRUE)

# ## CHOOSE 20 LANGUAGES RANDOMLY TO TEST THE CODE ##
# lang_sample <- df_values %>%
#   distinct(Language_ID) %>%
#   slice_sample(n = 10)
# # Filter the original dataset
# df_values <- df_values %>%
#   filter(Language_ID %in% lang_sample$Language_ID)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Run tests and write CSVs ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# NODF for r00 & c0
df_nodf_r00 <- nested_test_nodf_full(df_values, shuffling_type = "r00")
write.csv(df_nodf_r00, "nodf_results_all_lang_r00.csv", row.names = FALSE)

df_nodf_c0  <- nested_test_nodf_full(df_values, shuffling_type = "c0")
write.csv(df_nodf_c0, "nodf_results_all_lang_c0.csv", row.names = FALSE)



# Temperature for r00 & c0
df_temp_r00 <- nested_test_temp_full(df_values, shuffling_type = "r00")
write.csv(df_temp_r00, "temp_results_all_lang_r00.csv", row.names = FALSE)

df_temp_c0  <- nested_test_temp_full(df_values, shuffling_type = "c0")
write.csv(df_temp_c0, "temp_results_all_lang_c0.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. Distribution plot ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Function to plot the global Gaussian distributions
plot_distribution_global <- function(df_r00, df_c0,
                                     measure_label, x_label, x_lim) {
  # decide which column holds the values
  value_col <- if (tolower(measure_label) == "temperature") {
    "Value"
  } else {
    "NODF_global_Value"
  }
  # Extract simulated values
  sim_r00 <- df_r00 %>% filter(Type == "simulated") %>% pull(.data[[value_col]])
  sim_c0  <- df_c0  %>% filter(Type == "simulated") %>% pull(.data[[value_col]])
  # Compute mean and sd for each null model
  mean_r00 <- mean(sim_r00); sd_r00 <- sd(sim_r00)
  mean_c0  <- mean(sim_c0);  sd_c0  <- sd(sim_c0)
  # Extract observed (real) values
  real_r00 <- df_r00 %>% filter(Type == "real") %>% pull(.data[[value_col]])
  real_c0  <- df_c0  %>% filter(Type == "real") %>% pull(.data[[value_col]])
  # Create a grid for density curves
  x_seq <- seq(x_lim[1], x_lim[2], length.out = 300)
  df_density <- bind_rows(
    data.frame(x = x_seq,
               y = dnorm(x_seq, mean_r00, sd_r00),
               baseline = "r00"),
    data.frame(x = x_seq,
               y = dnorm(x_seq, mean_c0, sd_c0),
               baseline = "c0")
  )
  # Plot density curves and observed value lines
  ggplot(df_density, aes(x, y, fill = baseline, color = baseline)) +
    geom_area(alpha = 0.3, position = "identity") +
    geom_line(size = 1) +
    # add real‐value vertical lines
    geom_vline(
      xintercept = real_r00,
      linetype   = "dashed",
      size       = 1.1,
      color      = "red"
      # ) +
      # geom_vline(
      #   xintercept = real_c0,
      #   linetype   = "dashed",
      #   size       = 1.1,
      #   color      = "red"
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
    scale_x_continuous(limits = c(0,100)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
}


## Plot and save NODF distribution ----
dist_nodf_global <- plot_distribution_global(
  df_nodf_r00, df_nodf_c0,
  measure_label = "NODF",
  x_label       = "NODF",
  x_lim         = c(0, 100)
)
dist_nodf_global

# Save 
# ggsave("dist_nodf_global.png", dist_nodf_global,
#        width = 8, height = 6, bg = "white")

## Plot and save Temperature distribution ----
dist_temp_global <- plot_distribution_global(
  df_temp_r00, df_temp_c0,
  measure_label = "Temperature",
  x_label       = "Temperature",
  x_lim         = c(0, 100)
)
dist_temp_global
# Save 
# ggsave("dist_temp_global.png", dist_temp_global,
#        width = 8, height = 6, bg = "white")

