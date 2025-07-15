#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ==== Script for the analysis of the nestedness of phonological ====
# inventories across language families + initial analysis with SegBo 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This script contains 3 parts:
# I. Analysis of the nestedness of phonological family inventories
# II. Visualization of the inventory matrices themselves
# III. Analysis of borrowings with segbo

# Refer to the section header comments to identify the required input datasets 
# and the expected outputs of each section.

# REQUIRED INPUT FILES:
# 1. Phoible dataset (CLDF format):
#    - values.csv: contains phoneme inventories (phoneme occurrences)
#    - languages.csv: contains language metadata (family, Glottocode, etc.)
#
# 2. SegBo dataset:
#    - SegBo database - Phonemes.csv: contains information on borrowed phonemes 
#                                     across languages
#
# OPTIONAL (for loading precomputed results instead of running full simulations,
# CSVs files on GitHub):
#    - df_nodf_r00_new_value.csv, df_nodf_c0_new_value.csv   
#       => Simulated + real NODF values
#    - df_temp_r00.csv, df_temp_c0.csv  
#       => Simulated + real Temperature values

# ==== Libraries ====
library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(tidyr)
library(vegan)
library(permute)
library(lattice)
library(ggridges)
library(ggdist)
library(wesanderson)
library(patchwork)
library(ggtext)
library(tibble)
library(purrr)

theme_set(theme_bw())

# ==== Setting main parameters (alpha + nb of iterations) ====
ALPHA_ <- 0.005
N_ITER_ <- 1000

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# I. NESTEDNESS METRICS ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# In this section, we create the list of selected language families 
# (>= 4 languages) in order to build a phonological inventory binary matrix 
# for each family. 
# We then calculate the degree of nestedness for these matrices (families)
# using NODF and Temperature from the oecosimu package. 
# Finally, we use baselines and simulations (r00 and c0) to compare our results 
# with the vegan R package.

# This section can be run in one go and will compute both Temperature and NODF
# for the 70 matrices (families) included in the study.

# INPUT:
# - values.csv : phoneme values per language from Phoible (cldf file)
# - languages.csv : language metadata from Phoible (cldf file)
#
# OUTPUT:
# Nestedness results using NODF and Temperature
# - 4 CSV files containing simulated and real nestedness values:
#     * df_nodf_r00.csv, df_nodf_c0.csv
#     * df_temp_r00.csv, df_temp_c0.csv
# - 2 combined plots (r00 and c0) for NODF and Temperature:
#     * nodf_res_005.png, temp_res_005.png
# - Distribution of simulations plots for selected families (optional)
# - 2 CSV files obs_nodf.csv and obs_temp.csv with significance categories

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Useful functions ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ==== Function for nestedness testing ====
nested_test <- function(family_list, 
                        function_type = nestednodf, 
                        shuffling_type = 'r00',
                        save_ = TRUE,
                        n_iter = 2,
                        alt = "two.sided") {
  # Create the empty dataframe
  df_results <- data.frame(Family = character(), 
                           Measure = character(), 
                           Type = character(),
                           Value = numeric(),
                           p_value = numeric(),
                           stringsAsFactors = TRUE)
  
  for (fam_N in family_list) {
    ## Step 1: Select languages
    lgs <- df_languages[df_languages$Family_Name == fam_N, ]$ID
    ## Step 2: Select one inventory per language 
    selected_combinations <- df_values %>%
      filter(Language_ID %in% lgs) %>%
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
    ## Step 3: Retrieve matching rows
    df_selected <- df_values %>%
      dplyr::inner_join(selected_combinations, 
                        by = c("Language_ID", "Contribution_ID"))
    ## Step 4: One-hot encode phonemes
    df_one_hot <- df_selected %>%
      dplyr::select(Language_ID, Value) %>%
      dplyr::mutate(count = 1) %>%
      tidyr::pivot_wider(names_from = Value, values_from = count, 
                         values_fill = list(count = 0))
    ## Convert to matrix format
    df_one_hot_matrix <- as.matrix(df_one_hot[, -1])  # Remove Language_ID 
    rownames(df_one_hot_matrix) <- df_one_hot$Language_ID  # Set row names
    nrow_ <- nrow(df_one_hot_matrix)
    ncol_ <- ncol(df_one_hot_matrix)
    fill <- sum(df_one_hot_matrix) / (nrow_ * ncol_)
    ############ REMOVE FOR PROPER TESTING (after pre-reg, simulations) ########
    # df_one_hot_matrix <- matrix(rbinom(nrow_ * ncol_, 1, fill),
    #                             nrow = nrow_, ncol = ncol_,
    #                             dimnames = list(paste0("agent_",
    #                                                    LETTERS[1:nrow_]),
    #                                             paste0("item_", 1:ncol_)))
    ####
    print(fam_N)
    # Now, compute the simulations + the degree of nestedness using oecosimu 
    results <- oecosimu(df_one_hot_matrix, 
                        nestfun    = function_type, 
                        method     = shuffling_type,
                        parallel   = -1,
                        nsimul     = n_iter,
                        alternative = alt)
    # For NODF
    if (identical(function_type, nestednodf)) {
      # Select the third line of the vector (global NODF)
      row_idx     <- 3
      sim_global  <- results$oecosimu$simulated[row_idx, ] # Simulations nodf
      stat_global <- results$statistic$statistic[row_idx] # Real nodf
      pval_global <- results$oecosimu$pval[row_idx] # P value nodf
      n_sim       <- length(sim_global)
      # simulations rows
      new_rows <- data.frame(
        Family   = rep(fam_N, n_sim),
        Measure  = "NODF",
        Type     = "simulated",
        Value    = as.numeric(sim_global),
        p_value  = rep(as.numeric(pval_global), n_sim),
        n_langs  = rep(nrow(df_one_hot_matrix), n_sim)
      )
      # Real value row
      new_rows_real <- data.frame(
        Family   = fam_N,
        Measure  = "NODF",
        Type     = "real",
        Value    = as.numeric(stat_global),
        p_value  = as.numeric(pval_global),
        n_langs  = nrow(df_one_hot_matrix)
      )
    } else {
      sim <- results$oecosimu$simulated # Take all the simulations (simple vector)
      # For nested-temp
      # Simulations rows
      new_rows <- data.frame(
        Family = rep(fam_N, length(sim)),
        Measure = rep('Temperature', length(sim)),
        Type = rep('simulated', length(sim)),
        Value = as.numeric(sim),
        p_value = as.numeric(results$oecosimu$pval[1]), # Take the only p-value
        n_langs = nrow(df_one_hot_matrix)
      )
      # Real row
      new_rows_real <- data.frame(
        Family = fam_N,
        Measure = 'Temperature',
        Type = 'real',
        Value = as.numeric(results$statistic$statistic[1]), # Take the only value
        p_value = as.numeric(results$oecosimu$pval[1]), # Take the only p-value
        n_langs = nrow(df_one_hot_matrix)
      )
    }
    # Create the final dataframe
    df_results <- rbind(df_results, new_rows)
    df_results <- rbind(df_results, new_rows_real)
  }
  return(df_results)
}

## ==== Function to process simulated data (for both baselines) ====
process_simulated <- function(df_r00, df_c0) {
  df_sim_r00 <- df_r00 %>% filter(Type == "simulated") %>% 
    mutate(Baseline = "r00")
  df_sim_c0 <- df_c0 %>% filter(Type == "simulated") %>% 
    mutate(Baseline = "c0")
  # Combine simulated datasets into one
  df_sim <- bind_rows(df_sim_r00, df_sim_c0)
  # Compute quantiles for simulated data
  sim_summary <- df_sim %>% 
    group_by(Family, Baseline) %>% 
    summarise(lower = quantile(Value, 0.025),
              upper = quantile(Value, 0.975),
              .groups = "drop")
  return(list(df_sim = df_sim, sim_summary = sim_summary))
}

## ==== Function to process real data and create significance categories ====
process_real <- function(df_r00, df_c0, measure_label, alpha) {
  ## For r00
  df_obs_r00 <- df_r00 %>% group_by(Family) %>%
    # Calculate simulated mean per family
    mutate(Simulated_Mean = mean(Value)) %>%
    # Determine the significant side (Nested/Antinested)
    mutate(Significant_Side = case_when(
      Type != "real" ~ NA_character_,
      p_value > alpha ~ NA,
      Value > Simulated_Mean ~ 
        ifelse(measure_label == "NODF", "Nested", "Antinested"),
      Value < Simulated_Mean ~ 
        ifelse(measure_label == "NODF", "Antinested", "Nested"),
      TRUE ~ NA
    )) %>%
    filter(Type == "real") %>%
    # Assign a shape based on the direction of significance
    mutate(shape_type_r00 = case_when(
      Significant_Side == "Nested" ~ "triangle",   # Triangle for Nested
      Significant_Side == "Antinested" ~ "square", # Square for Antinested
      TRUE ~ "circle"                              # Circle if non significant
    )) %>%
    # Select the necessary columns
    select(Family, 
           n_langs,
           Value_r00 = Value,
           Simulated_Mean_r00 = Simulated_Mean,
           p_value_r00 = p_value,
           significant_r00 = significant,
           Significant_Side_r00 = Significant_Side,
           shape_type_r00)
  ## For c0
  df_obs_c0 <- df_c0 %>% group_by(Family) %>%
    # Calculate simulated mean per family
    mutate(Simulated_Mean = mean(Value)) %>%
    # Determine the significant side (Nested/Antinested)
    mutate(Significant_Side = case_when(
      Type != "real" ~ NA_character_,
      p_value >= alpha ~ NA,
      Value > Simulated_Mean ~ ifelse(measure_label == "NODF", 
                                      "Nested", "Antinested"),
      Value < Simulated_Mean ~ ifelse(measure_label == "NODF", 
                                      "Antinested", "Nested"),
      TRUE ~ NA
    )) %>%
    filter(Type == "real") %>%
    # Assign a shape based on the direction of significance
    mutate(shape_type_c0 = case_when(
      Significant_Side == "Nested" ~ "triangle",   # Triangle for Nested
      Significant_Side == "Antinested" ~ "square", # Square for Antinested
      TRUE ~ "circle"                              # Circle if non significant
    )) %>%
    # Select the necessary columns
    select(Family, 
           Value_c0 = Value, 
           Simulated_Mean_c0 = Simulated_Mean,
           p_value_c0 = p_value,
           significant_c0 = significant,
           Significant_Side_c0 = Significant_Side, 
           shape_type_c0)
  # Merge real data for both baselines in one dataframe
  df_obs_final <- inner_join(df_obs_r00, df_obs_c0, by = "Family") %>%
    # See if a family is significant for none, one or both baselines
    mutate(sig_cat = case_when(
      significant_r00 & significant_c0 ~ "both",
      xor(significant_r00, significant_c0) ~ "one",
      TRUE ~ "none"
    ))
  
  return(df_obs_final)
}

## ==== Function to plot combined (r00 and c0) observed data ====
plot_combined <- function(df_obs_final, sim_summary, x_label) {
  # Reorder families by the number of languages (ascending order)
  fam_order <- df_obs_final %>%
    distinct(Family, n_langs) %>%
    arrange(n_langs) %>%
    pull(Family)
  df_obs_final$Family <- factor(df_obs_final$Family, levels = fam_order)
  sim_summary$Family <- factor(sim_summary$Family, levels = fam_order) 
  # Bullet point for number of significant baselines per family
  # Significant for 0 baseline = black
  # Significant for 1 baseline = yellow
  # Significant for both baselines = green
  bullet_colors <- c("both" = "#009E73", "one" = "#FFA500", "none" = "black")
  axis_labels_df <- df_obs_final %>%
    distinct(Family, sig_cat) %>%
    mutate(label = paste0("<span style='color:black'>", Family, "</span> ",
                          "<span style='color:", bullet_colors[sig_cat],
                          "'>&#9679;</span>"))
  labels_vector <- setNames(axis_labels_df$label, axis_labels_df$Family)
  ## Plot for baseline r00
  p_r00 <- ggplot(df_obs_final, aes(x = Value_r00, y = Family)) +
    # Add segments representing simulated intervals for each family
    geom_segment(data = sim_summary %>% filter(Baseline == "r00"),
                 aes(x = lower, xend = upper, y = Family, yend = Family),
                 color = "#A9A9A9", linewidth = 0.8) +
    # Add observed points with fill and shape based on significance
    geom_point(aes(fill = shape_type_r00, shape = shape_type_r00),
               size = 3.5, color = "black") +
    # Set the different colors
    scale_fill_manual(name = "Significance type",
                      values = c("triangle" = "#64B5F6", # Nested
                                 "square" = "#FF6961", # Antinested
                                 "circle" = "grey"), # Non significant
                      labels = c("triangle" = "Nested",
                                 "square" = "Antinested",
                                 "circle" = "Not significant")) +
    # Set the different shapes
    scale_shape_manual(name = "Significance type",
                       values = c("triangle" = 24,
                                  "square" = 22,
                                  "circle" = 21),
                       labels = c("triangle" = "Nested",
                                  "square" = "Antinested",
                                  "circle" = "Not significant"))+
    # Set the axis and legend
    scale_y_discrete(labels = labels_vector) +
    labs(x = x_label, y = "", title = "r00") +
    scale_x_continuous(limits = c(0, 100)) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 9),
          axis.text.y = element_markdown(size = 10),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none")
  ## Plot for baseline c0
  p_c0 <- ggplot(df_obs_final, aes(x = Value_c0, y = Family)) +
    # Add segments representing simulated intervals for each family
    geom_segment(data = sim_summary %>% filter(Baseline == "c0"),
                 aes(x = lower, xend = upper, y = Family, yend = Family),
                 color = "#A9A9A9", linewidth = 0.8) +
    # Add observed points with fill and shape based on significance
    geom_point(aes(fill = shape_type_c0, shape = shape_type_c0),
               size = 3.5, color = "black") +
    # Set the different colors
    scale_fill_manual(name = "Significance type",
                      values = c("triangle" = "#64B5F6", # Nested
                                 "square" = "#FF6961", # Antinested
                                 "circle" = "grey"), # Non significant
                      labels = c("triangle" = "Nested",
                                 "square" = "Antinested",
                                 "circle" = "Not significant")) +
    # Set the different shapes
    scale_shape_manual(name = "Significance type",
                       values = c("triangle" = 24,
                                  "square" = 22,
                                  "circle" = 21),
                       labels = c("triangle" = "Nested",
                                  "square" = "Antinested",
                                  "circle" = "Not significant"))+
    # No legend for the y axis
    scale_y_discrete(labels = NULL) +
    # Set the axis and legend
    labs(x = x_label, title = "c0") +
    scale_x_continuous(limits = c(0, 100)) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 9),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none")
  # Combine (with patchwork)
  combined_plot <- (p_r00 + p_c0) +
    plot_layout(ncol = 2, guides = "collect") +
    plot_annotation(theme = theme(plot.title = element_text(hjust = 0.5))) &
    theme(legend.position = "bottom")
  return(combined_plot)
}



## ==== Function to plot the simulation distribution for a selected family ==== 
plot_distribution <- function(df_sim_r00, df_sim_c0, df_obs_r00, x_label,
                              selected_family, x_lim) {
  ########## CHANGE WHEN ENTIRE DATASET (NOT SELECTED FAMILY) #############
  df_sim_r00_entire_dataset <- df_sim_r00 %>% filter(Family == selected_family)
  df_sim_c0_entire_dataset  <- df_sim_c0 %>% filter(Family == selected_family)
  # Compute mean and standard deviation for each baseline
  mean_r00 <- mean(df_sim_r00_entire_dataset$Value, na.rm = TRUE)
  sd_r00   <- sd(df_sim_r00_entire_dataset$Value, na.rm = TRUE)
  mean_c0  <- mean(df_sim_c0_entire_dataset$Value, na.rm = TRUE)
  sd_c0    <- sd(df_sim_c0_entire_dataset$Value, na.rm = TRUE)
  # Extract real value
  real_r00 <- df_obs_r00 %>% filter(Family == selected_family) %>%
    pull(Value_r00)
  ####
  # Define x-axis range and density sequence
  x_seq <- seq(0, 100, length.out = 200)
  df_density_r00 <- data.frame(
    x = x_seq,
    y = dnorm(x_seq, mean = mean_r00, sd = sd_r00),
    baseline = "r00"
  )
  # Create data frames for the density curves
  df_density_c0 <- data.frame(
    x = x_seq,
    y = dnorm(x_seq, mean = mean_c0, sd = sd_c0),
    baseline = "c0"
  )
  # Create distribution plot
  p_gauss <- ggplot() +
    geom_area(data = df_density_r00, aes(x = x, y = y, fill = baseline),
              alpha = 0.3) +
    geom_area(data = df_density_c0, aes(x = x, y = y, fill = baseline),
              alpha = 0.3) +
    geom_line(data = df_density_r00, aes(x = x, y = y, color = baseline),
              linewidth = 1) +
    geom_line(data = df_density_c0, aes(x = x, y = y, color = baseline),
              linewidth = 1) +
    # r00 simulations in blue
    # c0 simulations in orange
    scale_fill_manual(values = c("r00" = "#56B", "c0" = "orange"),
                      guide = FALSE) +
    scale_color_manual(name = "Baseline",
                       values = c("r00" = "#56B", "c0" = "orange"),
                       labels = c("r00", "c0")) +
    labs(x = x_label, y = "Density",
         title = paste(selected_family)) +
    scale_x_continuous(limits = x_lim) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  # Add vertical line for the real value in red
  p_gauss <- p_gauss +
    geom_vline(xintercept = real_r00[1], color = "red", linetype = "dashed",
               linewidth = 1)
  return(p_gauss)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Loading and processing data ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ==== Read Phoible data ====
df_values <- read.csv("values.csv",
                      sep = ",",
                      header = TRUE)
df_languages <- read.csv("languages.csv",
                         sep = ",",
                         header = TRUE)

## ==== Generate list of families with at least 4 languages ====
f_n <- df_languages %>%
  group_by(Family_Name) %>%
  summarise(n_l = n()) %>%
  # Keep the families with more than 3 languages
  filter(n_l > 3)  %>%
  # Remove strange families
  filter(!Family_Name %in% c("", "Bookkeeping")) %>%
  pull(Family_Name)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Run Nestedness Tests ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ==== To import the datasets with the 1,000 simulations already in github ====
# df_temp_c0 <- read.csv("df_temp_c0.csv")
# df_temp_r00 <- read.csv("df_temp_r00.csv")
# df_nodf_c0 <- read.csv("df_nodf_c0_new_value.csv")
# df_nodf_r00 <- read.csv("df_nodf_r00_new_value.csv")

# ==== Run NODF tests ====
# 1) For r00
df_nodf_r00 <- nested_test(f_n, n_iter = N_ITER_, shuffling_type = "r00",
                           function_type = nestednodf)
# 2) For c0
df_nodf_c0  <- nested_test(f_n, n_iter = N_ITER_, shuffling_type = "c0",
                           function_type = nestednodf)

# 3) Save simulated datasets results NODF (simulations rows + real rows) 
write.csv(df_nodf_c0, "df_nodf_c0.csv", row.names = FALSE)
write.csv(df_nodf_r00, "df_nodf_r00.csv", row.names = FALSE)


## ==== Run Temperature tests ====
# 1) For r00
df_temp_r00 <- nested_test(f_n, n_iter = N_ITER_, shuffling_type = "r00",
                           function_type = nestedtemp)
# 2) For c0
df_temp_c0  <- nested_test(f_n, n_iter = N_ITER_, shuffling_type = "c0",
                           function_type = nestedtemp)

# 3) Save simulated datasets results temp (simulations rows + real rows) 
write.csv(df_temp_c0, "df_temp_c0.csv", row.names = FALSE)
write.csv(df_temp_r00, "df_temp_r00.csv", row.names = FALSE)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. Dataset Manipulation and Plotting ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ==== For NODF ====
# 1) Process simulated and observed NODF data (add significance)
df_nodf_r00_r <- df_nodf_r00 %>% mutate(significant = p_value <= ALPHA_)
df_nodf_c0_r  <- df_nodf_c0 %>% mutate(significant = p_value <= ALPHA_)
sim_nodf <- process_simulated(df_nodf_r00_r, df_nodf_c0_r)
obs_nodf <- process_real(df_nodf_r00_r, df_nodf_c0_r, "NODF", ALPHA_)

# 2) Create combined NODF plot 
plot_nodf <- plot_combined(obs_nodf, sim_nodf$sim_summary, "NODF")


# 3) Create NODF distribution plot for a selected family
dist_nodf <- plot_distribution(
  df_nodf_r00 %>% filter(Type == "simulated") %>% mutate(Baseline = "r00"),
  df_nodf_c0 %>% filter(Type == "simulated") %>% mutate(Baseline = "c0"),
  obs_nodf, "NODF", "Algic", c(20, 65)
)


## ==== For Temperature ====
# 1) Process simulated and observed Temperature data (add significance)
df_temp_r00_r <- df_temp_r00 %>% mutate(significant = p_value <= ALPHA_)
df_temp_c0_r  <- df_temp_c0 %>% mutate(significant = p_value <= ALPHA_)
sim_temp <- process_simulated(df_temp_r00_r, df_temp_c0_r)
obs_temp <- process_real(df_temp_r00_r, df_temp_c0_r, "Temperature", ALPHA_)

# 2) Create combined Temperature plot
plot_temp <- plot_combined(obs_temp, sim_temp$sim_summary, "Temperature")
plot_temp 

# 3) Create Temperature distribution plot for a selected family
dist_temp <- plot_distribution(
  df_temp_r00 %>% filter(Type == "simulated") %>% mutate(Baseline = "r00"),
  df_temp_c0 %>% filter(Type == "simulated") %>% mutate(Baseline = "c0"),
  obs_temp, "Temperature", "Algic", c(20,65)
)


# ==== Save final dataset results for temp and NODF ====
write.csv(obs_nodf, "obs_nodf.csv", row.names = FALSE)
write.csv(obs_temp, "obs_temp.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. Display the Plots ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ==== Display the plots ====
plot_nodf 
plot_temp

## ==== Save combined plots for NODF and Temperature ====
ggsave('nodf_res_005.png',
       plot=plot_nodf,
       width=14,
       height=10,
       scale=0.9,
       bg='white')
ggsave('temp_res_005.png',
       plot=plot_temp,
       width=14,
       height=10,
       scale=0.9,
       bg='white')




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  II. MATRIX ANALYSIS AND VISUALISATION ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This section displays the phonological inventory matrices (real or simulated).
# Some of these functions / bits of code are not essential for the final analysis
# and can be deleted 

# INPUT (already loaded in section I):
# - values.csv : phoneme values per language from Phoible (cldf file)
# - languages.csv : language metadata from Phoible (cldf file)
#
# OUTPUT:
# - Sorted and visualized inventory matrices (real and simulated)
# - PNG figures for:
#     * Big families' matrices (panel_big_families_matrices.png)
#     * Real + simulated matrix panels (panel_Mande_r00.png / panel_Mande_c0.png)
# - Side-by-side visual comparison (original vs sorted matrices)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Function Definitions ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ==== Function to sort a matrix by rows and columns sums ====
sort_matrix <- function(M, iterations = 100) {
  # Perform iterative sorting to stabilize the order
  for (i in seq_len(iterations)) {
    # Sort columns in decreasing order of their sums (most 1s first)
    if (ncol(M) > 1) {
      M <- M[, order(colSums(M), decreasing = TRUE)]
    }
    # Sort rows in decreasing order of their sums (most 1s first)
    if (nrow(M) > 1) {
      M <- M[order(rowSums(M), decreasing = TRUE), ]
    }
  }
  return(M)
}

## ==== Function to plot a matrix ====
plot_matrix <- function(mat, title = "") {
  # Convert the matrix to long format for ggplot
  df <- melt(mat)
  # Rename the columns for clarity
  colnames(df) <- c("Language", "Phoneme", "Value")
  # Reverse the order of languages so they display top-to-bottom in the plot
  df$Language <- factor(df$Language, levels = rev(rownames(mat)))
  # Create the plot
  p <- ggplot(df, aes(x = Phoneme, y = Language, fill = factor(Value))) +
    # Draw tiles for each matrix cell
    geom_tile(color= "white") +
    # colors: 0 = beige, 1 = black
    scale_fill_manual(values = c("0" = "cornsilk", "1" = "black"),
                      guide = "none") +
    labs(x = "", y = "", title = title) +
    theme_minimal() +
    theme(axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0.5))
  return(p)
}


## ==== Function to extract the binary matrix for a family ====
get_family_matrix <- function(family_name, 
                              values = df_values, 
                              languages = df_languages) {
  ## Step 1: Select languages (the same as for the nestedness analysis in I.)
  lgs <- languages$ID[languages$Family_Name == family_name]
  ## Step 2: Select one inventory per language 
  selected_contributions <- values %>%
    filter(Language_ID %in% lgs) %>%
    group_by(Language_ID, Contribution_ID) %>%
    summarise(
      phoneme_list = list(unique(Value)),
      n_phonemes = n_distinct(Value),
      .groups = 'drop'
    ) %>%
    group_by(Language_ID) %>%
    group_modify(~ {
      if (nrow(.x) == 1) {
        return(.x)
      } else if (nrow(.x) == 2) {
        min_count <- min(.x$n_phonemes)
        candidates <- .x %>% filter(n_phonemes == min_count)
        return(candidates %>% slice_sample(n = 1))
      } else {
        overlaps <- sapply(1:nrow(.x), function(i) {
          sum(sapply(1:nrow(.x), function(j) {
            if (i != j) length(intersect(.x$phoneme_list[[i]], 
                                         .x$phoneme_list[[j]]))
            else 0
          }))
        })
        candidates <- .x[overlaps == max(overlaps), ]
        min_count <- min(candidates$n_phonemes)
        candidates <- candidates %>% filter(n_phonemes == min_count)
        return(candidates %>% slice_sample(n = 1))
      }
    }) %>%
    ungroup() %>%
    select(Language_ID, Contribution_ID)
  ## Step 3: Retrieve matching rows
  df_selected <- values %>% 
    inner_join(selected_contributions, by = c("Language_ID", 
                                              "Contribution_ID"))
  ## Step 4: One-hot encode phonemes
  df_one_hot <- df_selected %>% 
    select(Language_ID, Value) %>% 
    mutate(count = 1) %>% 
    pivot_wider(names_from = Value, values_from = count, 
                values_fill = list(count = 0))
  ## Convert to matrix format
  binary_matrix <- as.matrix(df_one_hot[, -1])
  rownames(binary_matrix) <- df_one_hot$Language_ID
  return(binary_matrix)
}

## ==== Function to simulate matrices using r00 ====
simulate_r00 <- function(real_mat, n_sim = 3, sort_iter = 1000) {
  # Calculate the fill probability
  fill_prob <- sum(real_mat) / (nrow(real_mat) * ncol(real_mat))
  # Initialize a list to store the plot of each simulated matrix
  sim_plots <- vector("list", n_sim)
  # Loop over the number of simulations (n_sim)
  for (i in 1:n_sim) {
    # For each cell, randomly generate a 0 or a 1,
    # where the chance of getting a 1 is given by fill_prob
    sim_mat <- matrix(rbinom(nrow(real_mat) * ncol(real_mat), 
                             size = 1, 
                             prob = fill_prob),
                      nrow = nrow(real_mat), ncol = ncol(real_mat),
                      dimnames = list(rownames(real_mat), colnames(real_mat)))
    # Sort the simulated matrix
    sim_sorted <- sort_matrix(sim_mat, iterations = sort_iter)
    # Create a plot of the sorted simulated matrix
    sim_plots[[i]] <- plot_matrix(sim_sorted, paste("Simulated matrix (r00)", i))
  }
  # Return the list of simulated matrix plots
  return(sim_plots)
}


## ==== Function to simulate matrices using c0 ====
simulate_c0 <- function(real_mat, n_sim = 3, sort_iter = 1000) {
  # calculate the probability for each column (fraction of 1s in that column)
  col_probs <- colSums(real_mat) / nrow(real_mat)
  # Create an empty list to store the plots of the simulated matrices
  sim_plots <- vector("list", n_sim)
  # Loop over the number of simulations
  for (i in 1:n_sim) {
    # create an empty matrix with the same dimensions as the observed matrix
    sim_mat <- matrix(NA, nrow = nrow(real_mat), 
                      ncol = ncol(real_mat),
                      dimnames = list(rownames(real_mat), 
                                      colnames(real_mat)))
    # For each column generate 0 or 1 randomly
    # the probability of getting a 1 in that column is given by col_probs[j]
    for (j in seq_along(col_probs)) {
      sim_mat[, j] <- rbinom(n = nrow(real_mat), size = 1, prob = col_probs[j])
    }
    # sort the simulated matrix
    sim_sorted <- sort_matrix(sim_mat, iterations = sort_iter)
    # Create a plot of the sorted simulated matrix
    sim_plots[[i]] <- plot_matrix(sim_sorted, paste("Simulated matrix (c0)", i))
  }
  # Return the list of simulated matrix plots
  return(sim_plots)
}

## ==== Function to arrange and display a panel of plots (real + simulated) ====
plot_panel <- function(real_plot, sim_plots, ncol = 2) {
  # Combine real plot with simulated plots
  all_plots <- c(list(real_plot), sim_plots)
  # Arrange plots in a grid with specified columns
  panel <- arrangeGrob(grobs = all_plots, ncol = ncol)
  # Display the panel
  grid.arrange(panel)
  return(panel)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Panel for big families matrices ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ==== Select the biggest families ====
selected_big_families <- families_class$Family_Name[1:6]
print(selected_big_families)

# ==== For each selected family, extract, sort, and plot the binary matrix ====
big_families_plots <- lapply(selected_big_families, function(fam) {
  bin_mat <- get_family_matrix(fam)
  sorted_mat <- sort_matrix(bin_mat, iterations = 1000)
  plot_matrix(sorted_mat, fam)
})

# ==== Create a panel ====
panel_big_families <- arrangeGrob(grobs = big_families_plots, ncol = 3)
grid.arrange(panel_big_families)

# ==== Save the panel ====
ggsave(filename = "panel_big_families_matrices.png",
       plot = panel_big_families,
       width = 20,
       height = 12, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Display a single family matrix (unsorted and sorted) ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ==== Select the family (here "Uralic") ====
uralic_matrix <- get_family_matrix("Uralic")

# ==== Plot the original (unsorted) inventory matrix ====
p_original_uralic <- plot_matrix(uralic_matrix, "Original Uralic Matrix")
# print(p_original_uralic)

# ==== Plot the sorted matrix ====
sorted_uralic_matrix <- sort_matrix(uralic_matrix, iterations = 1000)
p_sorted_uralic <- plot_matrix(sorted_uralic_matrix, "Sorted Uralic Matrix")
# print(p_sorted_uralic)

# ==== To display both plots side by side ====
p_original_uralic + p_sorted_uralic


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. Panel of the simulated matrices (r00 and c0) ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ==== Extract and sort the real matrix ====
real_mat_big <- get_family_matrix("Mande")
real_sorted_big <- sort_matrix(real_mat_big, iterations = 10000)
p_real_big <- plot_matrix(real_sorted_big, "Mande")

## ==== Simulate matrices using the r00 model ====
set.seed(123)
sim_plots_big_r00 <- simulate_r00(real_mat_big,
                                  n_sim = 3,
                                  sort_iter = 1000)
panel_big_r00 <- plot_panel(p_real_big,
                            sim_plots_big_r00,
                            ncol = 2)

# ==== To save the panel ====
ggsave(filename = "panel_Mande_r00.png",
       plot = panel_big_r00,
       width = 12, height = 10,
       dpi = 300)

## ==== Simulate matrices using the c0 model ====
set.seed(123)
sim_plots_big_c0 <- simulate_c0(real_mat_big,
                                n_sim = 3,
                                sort_iter = 1000)
p_real_big_c0 <- plot_matrix(real_sorted_big, "Mande")
panel_big_c0 <- plot_panel(p_real_big_c0,
                           sim_plots_big_c0,
                           ncol = 2)

# ==== To save the panel ====
ggsave(filename = "panel_Mande_c0.png",
       plot = panel_big_c0,
       width = 12,
       height = 10,
       dpi = 300)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# III. OVERLAP SEGBO ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This section performs the overlap between the dataset used for the nestedness 
# analysis in section I and the SegBo dataset, allowing us to identify borrowed 
# phonemes within phonological inventories.

# This analysis is not complete, and this code mainly serves to determine 
# the number of borrowings in the overlap, as well as to gather additional 
# information, without deeply analyzing the contribution (positive or negative) 
# of borrowings to nestedness.

# INPUT:
# - SegBo dataset: "SegBo database - Phonemes.csv"
# - Phoible values and languages (values.csv and languages.csv, 
#   already loaded in section I)
#
# OUTPUT:
# - Visualization of borrowed phonemes in inventories:
#     * matrix_borrowings_Austronesian.png
# - Overlap dataset of borrowed phonemes across major families
# - Summary tables of borrowing counts per family
# - Analysis of marginal phonemes in Phoible vs SegBo presence

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Load SegBo dataset ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ==== Import segbo dataset ====
df_segbo <- read.csv(
  "SegBo database - Phonemes.csv",
  sep  = ",",
  header   = TRUE,
  stringsAsFactors  = FALSE
) %>%
  rename(
    language = BorrowingLanguageGlottocode,
    feature  = BorrowedSound,
    only_loanwords  = OnlyInLoanwords,
    result  = Result,
    new_distinction = NewDistinction,
    comments  = PhonemeComments,
    verified   = Verified
  ) %>%
  mutate(borrowed = 1L)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Functions ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ==== Function to plot a matrix with its borrowed phonemes ====
plot_family_borrowed <- function(family_name, df_segbo, sort_iter = 1000) {
  # Extract the matrix and sort it
  mat        <- get_family_matrix(family_name)
  mat_sorted <- sort_matrix(mat, iterations = sort_iter)
  # Melt the sorted matrix into long format
  df <- melt(mat_sorted,
             varnames  = c("language", "feature"),
             value.name = "present") %>%
    # Join with segbo data to detect borrowings
    left_join(df_segbo, by = c("language", "feature")) %>%
    # categorical status for plotting
    mutate(
      status = factor(
        case_when(
          present == 0 ~ "absent", 
          # absent if phoneme not in inventory
          present == 1 & borrowed == 1 ~ "borrowed",
          # borrowed if in inventory and flagged by SegBo
          TRUE ~ "present"
          # present if in inventory but not marked as borrowed
        ),
        levels = c("absent", "present", "borrowed")
      ),
      # align axis order with sorted matrix
      language = factor(language, levels = rev(rownames(mat_sorted))),
      feature  = factor(feature,  levels = colnames(mat_sorted))
    )
  # Create the plot and add colors
  ggplot(df, aes(x = feature, y = language, fill = status)) +
    geom_tile(color = "white") +
    scale_fill_manual(
      values = c("absent" = "cornsilk",
                 "present" = "grey",
                 "borrowed"= "red"), # borrowed phonemes in red
      guide = "none"
    ) +
    # labels and theme
    labs(
      title = paste0(family_name, " – Inventory & SegBo Borrowings"),
      x = NULL, y = NULL
    ) +
    theme_minimal() +
    theme(
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
}

## ==== Function for various datasets about borrowings for a given family ====
get_family_data <- function(family_name, sort_iter = 1000) {
  # 1) Extract the matrix and sort it
  mat        <- get_family_matrix(family_name)
  mat_sorted <- sort_matrix(mat, iterations = sort_iter)
  # 2) Convert the sorted matrix to a data frame 
  df_matrix <- as.data.frame(mat_sorted) %>%
    rownames_to_column(var = "language")
  # 3) Pivot to long and keep only present phonemes
  df_long <- df_matrix %>%
    pivot_longer(
      cols = -language,
      names_to  = "feature",
      values_to = "present"
    ) %>%
    filter(present == 1) %>%
    select(-present)
  # 4) Join with SegBo to get full borrowing info
  merge_segbo <- df_long %>%
    inner_join(df_segbo, by = c("language", "feature")) %>%
    mutate(borrowed_flag = 1L)
  # 5) Extract unique borrowed (language, feature) pairs
  segbo_inter <- merge_segbo %>%
    filter(borrowed_flag == 1L) %>%
    distinct(language, feature)
  # 6) Return a list of all intermediate tables
  list(
    df_matrix   = df_matrix,
    df_long     = df_long,
    merge_segbo = merge_segbo,
    segbo_inter = segbo_inter
  )
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Get segbo data overlap for a family ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ==== Select a family and get its matrix ====
res_austronesian   <- get_family_data("Austronesian")
# 1) Sorted binary presence matrix
austronesian_df       <- res_austronesian$df_matrix
# 2) long-format (present==1)
austronesian_long_present<- res_austronesian$df_long
# 3) Full join with SegBo
austronesian_merge_segbo <- res_austronesian$merge_segbo
# 4) keep only (language, features)
austronesian_segbo_inter <- res_austronesian$segbo_inter

## ==== Plot the matrix with its borrowed phonemes ====
plot_austronesian <- plot_family_borrowed("Austronesian", df_segbo)
print(plot_austronesian)

## ==== To save the plot ====
ggsave(filename = "matrix_borrowings_Austronesian.png",
       plot = plot_austronesian,
       width = 10,
       height = 8,
       dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. Create the full dataset overlap Segbo / Phoible ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ==== (Re)load Phoible data if necessary ====
df_values <- read.csv("values.csv",
                      sep = ",",
                      header = TRUE)
df_languages <- read.csv("languages.csv",
                         sep = ",",
                         header = TRUE)

## ==== list of families that have at least 10 languages ====
## (otherwise not enough statistical power to study borrowings)
fam_list <- df_languages %>%
  filter(Family_Name != "", Family_Name != "Bookkeeping") %>%
  count(Family_Name) %>%
  filter(n > 9) %>%
  pull(Family_Name)

## ==== Create the overlap dataset ====
all_segbo_inter <- map_dfr(fam_list, function(fam) {
  # Retrieve and sort the binary phoneme matrix for this family
  mat <- get_family_matrix(fam, values = df_values, languages = df_languages)
  # long format
  df_long <- as.data.frame(mat) %>%
    rownames_to_column("language") %>%
    pivot_longer(cols = -language,
                 names_to  = "feature",
                 values_to = "present") %>%
    filter(present == 1) %>%
    select(language, feature)
  # Join with SegBo 
  df_long %>%
    inner_join(df_segbo, by = c("language", "feature")) %>%
    mutate(family = fam)
})

## ==== Compute number of languages per family ====
family_sizes <- df_languages %>%
  filter(Family_Name != "", Family_Name != "Bookkeeping") %>%
  count(Family_Name, name = "n_languages")

## ==== Summarize borrowings by family sorting by family size ====
# count borrowings per family
borrow_counts <- all_segbo_inter %>%
  count(family, name = "n_borrowings")
# select families with ≥10 languages
big_families <- family_sizes %>%
  filter(n_languages >= 10) %>%
  rename(family = Family_Name)
# include all families, set missing borrowings to 0
fam_overlap <- big_families %>%
  left_join(borrow_counts, by = "family") %>%
  mutate(n_borrowings = replace_na(n_borrowings, 0L)) %>%
  arrange(desc(n_languages))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. Marginal phonemes in Phoible vs segbo presence ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# We want to study the phonemes marked as marginal in the phoible dataset 
# and whether they are present in the segbo dataset.

# turn the marginal of phoible column into TRUE/FALSE
df_values$Marginal <- as.logical(df_values$Marginal)

## ==== Extract all marginal phonemes from Phoible ====
marginal_phoible_df <- df_values %>%
  filter(Marginal) %>%
  transmute(language = Language_ID,
            feature  = Value) 

## Those that *are* in SegBo
marginal_in_segbo <- marginal_phoible_df %>%
  inner_join(df_segbo, by = c("language", "feature"))

## Those that *aren’t* in SegBo
marginal_not_in_segbo <- marginal_phoible_df %>%
  anti_join(df_segbo, by = c("language", "feature"))

## Counts
total_marginal      <- nrow(marginal_phoible_df)
in_segbo_count      <- nrow(marginal_in_segbo)
not_in_segbo_count  <- nrow(marginal_not_in_segbo)

cat("Total marginal phonemes in Phoible: ", total_marginal,      "\n",
    "Found in SegBo:                        ", in_segbo_count,      "\n",
    "Not found in SegBo (dropped rows):    ", not_in_segbo_count,  "\n")

## See how phonemes called "marginal" in Phoible are annotated 
## in SegBo’s OnlyInLoanwords column 
## (if they truly are marginal, the value should be "yes")
marginal_loanwords_counts <- marginal_in_segbo %>%
  count(only_loanwords, name = "n") %>%
  arrange(desc(n))

