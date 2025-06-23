# Libraries
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

# Setting main parameters
ALPHA_ <- 0.005
N_ITER_ <- 1000

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# I. NESTEDNESS METRICS ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Useful functions ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Function for nestedness testing
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
    ############ REMOVE FOR PROPER TESTING (after pre-reg) ###############
    # df_one_hot_matrix <- matrix(rbinom(nrow_ * ncol_, 1, fill),
    #                             nrow = nrow_, ncol = ncol_,
    #                             dimnames = list(paste0("agent_",
    #                                                    LETTERS[1:nrow_]),
    #                                             paste0("item_", 1:ncol_)))
    ####
    print(fam_N)
    results <- oecosimu(df_one_hot_matrix, 
                        nestfun = function_type, 
                        method = shuffling_type,
                        parallel = -1,
                        nsimul = n_iter,
                        alternative = alt)
    sim <- results$oecosimu$simulated
    if (identical(function_type, nestednodf)) {
      # Select the third line (global NODF)
      row_idx     <- 3
      sim_global  <- results$oecosimu$simulated[row_idx, ]
      stat_global <- results$statistic$statistic[row_idx]
      pval_global <- results$oecosimu$pval[row_idx]
      n_sim       <- length(sim_global)
      # simulations
      new_rows <- data.frame(
        Family   = rep(fam_N, n_sim),
        Measure  = "NODF",
        Type     = "simulated",
        Value    = as.numeric(sim_global),
        p_value  = rep(as.numeric(pval_global), n_sim),
        n_langs  = rep(nrow(df_one_hot_matrix), n_sim)
      )
      # Real Value
      new_rows_real <- data.frame(
        Family   = fam_N,
        Measure  = "NODF",
        Type     = "real",
        Value    = as.numeric(stat_global),
        p_value  = as.numeric(pval_global),
        n_langs  = nrow(df_one_hot_matrix)
      )
    } else {
      sim <- results$oecosimu$simulated
      # For nested-temp
      new_rows <- data.frame(
        Family = rep(fam_N, length(sim)),
        Measure = rep('Temperature', length(sim)),
        Type = rep('simulated', length(sim)),
        Value = as.numeric(sim),
        p_value = as.numeric(results$oecosimu$pval[1]),
        n_langs = nrow(df_one_hot_matrix)
      )
      new_rows_real <- data.frame(
        Family = fam_N,
        Measure = 'Temperature',
        Type = 'real',
        Value = as.numeric(results$statistic$statistic[1]),
        p_value = as.numeric(results$oecosimu$pval[1]),
        n_langs = nrow(df_one_hot_matrix)
      )
    }
    df_results <- rbind(df_results, new_rows)
    df_results <- rbind(df_results, new_rows_real)
  }
  return(df_results)
}
## Function to process simulated data (for both baselines)
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

## Function to process real observed data and create significance categories
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

## Function to plot combined (r00 and c0) observed data
plot_combined <- function(df_obs_final, sim_summary, x_label) {
  # Reorder families by the number of languages (ascending order)
  fam_order <- df_obs_final %>%
    distinct(Family, n_langs) %>%
    arrange(n_langs) %>%
    pull(Family)
  df_obs_final$Family <- factor(df_obs_final$Family, levels = fam_order)
  sim_summary$Family <- factor(sim_summary$Family, levels = fam_order) 
  # Bullet point for number of significant baselines per family
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
                      values = c("triangle" = "#64B5F6",
                                 "square" = "#FF6961",
                                 "circle" = "grey"),
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
                      values = c("triangle" = "#64B5F6",
                                 "square" = "#FF6961",
                                 "circle" = "grey"),
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
  # Combine with patchwork
  combined_plot <- (p_r00 + p_c0) +
    plot_layout(ncol = 2, guides = "collect") +
    plot_annotation(theme = theme(plot.title = element_text(hjust = 0.5))) &
    theme(legend.position = "bottom")
  return(combined_plot)
}


## Function to plot the Gaussian distribution for the entire dataset
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
  # Add vertical line for the real value
  p_gauss <- p_gauss +
    geom_vline(xintercept = real_r00[1], color = "red", linetype = "dashed",
               linewidth = 1)
  return(p_gauss)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Loading and processing data ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Read Phoible data
df_values <- read.csv("data/cldf-datasets-phoible-f36deac/cldf/values.csv",
                      sep = ",",
                      header = TRUE)
df_languages <- read.csv("data/cldf-datasets-phoible-f36deac/cldf/languages.csv",
                         sep = ",",
                         header = TRUE)


## Generate list of families
f_n <- df_languages %>%
  group_by(Family_Name) %>%
  summarise(n_l = n()) %>%
  filter(n_l > 3)  %>%
  filter(!Family_Name %in% c("", "Bookkeeping")) %>%
  pull(Family_Name)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. VOWEL: Run Nestedness Tests ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Select only the vowels
df_params <- read.csv("data/parameters.csv",
                      header = TRUE, sep = ",")
df_vowel <- df_params %>% 
  filter(SegmentClass == "vowel") %>% 
  select(ID)

## Keep only the vowel in df_values
df_values <- df_values %>%
  inner_join(df_vowel, by = c("Parameter_ID" = "ID")) %>%
  select(Language_ID, Contribution_ID, Value)



# Run NODF tests
# For r00
df_nodf_vow_r00 <- nested_test(f_n, n_iter = N_ITER_, shuffling_type = "r00",
                           function_type = nestednodf)
# For c0
df_nodf_vow_c0  <- nested_test(f_n, n_iter = N_ITER_, shuffling_type = "c0",
                           function_type = nestednodf)

## Run Temperature tests

# # Remove "Koman" and "Maningrida" for temp metric (full matrices)
# f_n <- f_n[!f_n %in% c("Koman", "Maningrida")]
# 
# # For r00
# df_temp_vow_r00 <- nested_test(f_n, n_iter = N_ITER_, shuffling_type = "r00",
#                            function_type = nestedtemp)
# # For c0
# df_temp_vow_c0  <- nested_test(f_n, n_iter = N_ITER_, shuffling_type = "c0",
#                            function_type = nestedtemp)

# df_temp_c0 <- read.csv("data/df_temp_c0.csv")
# df_temp_vow_r00 <- read.csv("Dataset consonne vowel/Vowel/df_temp_r00_vow_1000.csv")
# df_nodf_vow_c0 <- read.csv("Dataset consonne vowel/Vowel/df_nodf_c0_vow_1000.csv")
# df_nodf_vow_r00 <- read.csv("Dataset consonne vowel/Vowel/df_nodf_r00_vow_1000.csv")


# # Save simulated datasets results for temp and NODF
# write.csv(df_temp_vow_c0, "df_temp_c0_vow_1000.csv", row.names = FALSE)
# write.csv(df_temp_vow_r00, "df_temp_r00_vow_1000.csv", row.names = FALSE)
write.csv(df_nodf_vow_c0, "df_nodf_c0_vow_1000_new_nodf.csv", row.names = FALSE)
write.csv(df_nodf_vow_r00, "df_nodf_r00_vow_1000_new_nodf.csv", row.names = FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. VOWEL: Dataset Manipulation and Plotting ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## For NODF
# Process simulated and observed NODF data
df_nodf_vow_r00_r <- df_nodf_vow_r00 %>% mutate(significant = p_value <= ALPHA_)
df_nodf_vow_c0_r  <- df_nodf_vow_c0 %>% mutate(significant = p_value <= ALPHA_)
sim_nodf_vow <- process_simulated(df_nodf_vow_r00_r, df_nodf_vow_c0_r)
obs_nodf_vow <- process_real(df_nodf_vow_r00_r, df_nodf_vow_c0_r, "NODF", ALPHA_)
# Create combined NODF plot
plot_nodf_vow <- plot_combined(obs_nodf_vow, sim_nodf_vow$sim_summary, "NODF")
# Create NODF distribution plot
dist_nodf <- plot_distribution(
  df_nodf_vow_r00 %>% filter(Type == "simulated") %>% mutate(Baseline = "r00"),
  df_nodf_vow_c0 %>% filter(Type == "simulated") %>% mutate(Baseline = "c0"),
  obs_nodf_vow, "NODF", "Algic", c(20, 65)
)

# ## For Temperature
# # Process simulated and observed Temperature data
# df_temp_vow_r00_r <- df_temp_vow_r00 %>% mutate(significant = p_value <= ALPHA_)
# df_temp_vow_c0_r  <- df_temp_vow_c0 %>% mutate(significant = p_value <= ALPHA_)
# sim_temp_vow <- process_simulated(df_temp_vow_r00_r, df_temp_vow_c0_r)
# obs_temp_vow <- process_real(df_temp_vow_r00_r, df_temp_vow_c0_r, "Temperature", ALPHA_)
# # Create combined Temperature plot
# plot_temp_vow <- plot_combined(obs_temp_vow, sim_temp_vow$sim_summary, "Temperature")
# # Create Temperature distribution plot
# dist_temp <- plot_distribution(
#   df_temp_vow_r00 %>% filter(Type == "simulated") %>% mutate(Baseline = "r00"),
#   df_temp_vow_c0 %>% filter(Type == "simulated") %>% mutate(Baseline = "c0"),
#   obs_temp_vow, "Temperature", "Algic", c(20,65)
# )


# # Save final dataset results for temp and NODF
# write.csv(obs_nodf_vow, "obs_nodf_vow_1000.csv", row.names = FALSE) 
# write.csv(obs_temp_vow, "obs_temp_vow_1000.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. VOWEL: Display the Plots ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Print combined plots for NODF and Temperature

ggsave('nodf_res_1000_vow_new_nodf.png',
       plot=plot_nodf_vow,
       width=14,
       height=10,
       scale=0.9,
       bg='white')
# ggsave('temp_res_vow_1000_new_nodf.png',
#        plot=plot_temp_vow,
#        width=14,
#        height=10,
#        scale=0.9,
#        bg='white')


# plot_temp_vow
plot_nodf_vow



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6. CONSONANTS: Run Nestedness Tests ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Read Phoible data
df_values <- read.csv("data/cldf-datasets-phoible-f36deac/cldf/values.csv",
                      sep = ",",
                      header = TRUE)

## Select only the consonants
df_params <- read.csv("data/parameters.csv",
                      header = TRUE, sep = ",")
df_consonant <- df_params %>% 
  filter(SegmentClass == "consonant") %>% 
  select(ID)

## Keep only the consonants in df_values
df_values <- df_values %>%
  inner_join(df_consonant, by = c("Parameter_ID" = "ID")) %>%
  select(Language_ID, Contribution_ID, Value)



# Run NODF tests
# For r00
df_nodf_cons_r00 <- nested_test(f_n, n_iter = N_ITER_, shuffling_type = "r00",
                               function_type = nestednodf)
# For c0
df_nodf_cons_c0  <- nested_test(f_n, n_iter = N_ITER_, shuffling_type = "c0",
                               function_type = nestednodf)

## Run Temperature tests

# # For r00
# df_temp_cons_r00 <- nested_test(f_n, n_iter = N_ITER_, shuffling_type = "r00",
#                                function_type = nestedtemp)
# # For c0
# df_temp_cons_c0  <- nested_test(f_n, n_iter = N_ITER_, shuffling_type = "c0",
#                                function_type = nestedtemp)



# Save simulated datasets results for temp and NODF
# write.csv(df_temp_cons_c0, "df_temp_c0_cons_1000.csv", row.names = FALSE)
# write.csv(df_temp_cons_r00, "df_temp_r00_cons_1000.csv", row.names = FALSE)
write.csv(df_nodf_cons_c0, "df_nodf_c0_cons_1000_new_nodf.csv", row.names = FALSE)
write.csv(df_nodf_cons_r00, "df_nodf_r00_cons_1000_new_nodf.csv", row.names = FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 7. CONSONANTS: Dataset Manipulation and Plotting ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## For NODF
# Process simulated and observed NODF data
df_nodf_cons_r00_r <- df_nodf_cons_r00 %>% mutate(significant = p_value <= ALPHA_)
df_nodf_cons_c0_r  <- df_nodf_cons_c0 %>% mutate(significant = p_value <= ALPHA_)
sim_nodf_cons <- process_simulated(df_nodf_cons_r00_r, df_nodf_cons_c0_r)
obs_nodf_cons <- process_real(df_nodf_cons_r00_r, df_nodf_cons_c0_r, "NODF", ALPHA_)
# Create combined NODF plot
plot_nodf_cons <- plot_combined(obs_nodf_cons, sim_nodf_cons$sim_summary, "NODF")
# Create NODF distribution plot
dist_nodf <- plot_distribution(
  df_nodf_cons_r00 %>% filter(Type == "simulated") %>% mutate(Baseline = "r00"),
  df_nodf_cons_c0 %>% filter(Type == "simulated") %>% mutate(Baseline = "c0"),
  obs_nodf_cons, "NODF", "Algic", c(20, 65)
)

# ## For Temperature
# # Process simulated and observed Temperature data
# df_temp_cons_r00_r <- df_temp_cons_r00 %>% mutate(significant = p_value <= ALPHA_)
# df_temp_cons_c0_r  <- df_temp_c0 %>% mutate(significant = p_value <= ALPHA_)
# sim_temp_cons <- process_simulated(df_temp_cons_r00_r, df_temp_cons_c0_r)
# obs_temp_cons <- process_real(df_temp_cons_r00_r, df_temp_cons_c0_r, "Temperature", ALPHA_)
# # Create combined Temperature plot
# plot_temp_cons <- plot_combined(obs_temp_cons, sim_temp_cons$sim_summary, "Temperature")
# # Create Temperature distribution plot
# dist_temp <- plot_distribution(
#   df_temp_cons_r00 %>% filter(Type == "simulated") %>% mutate(Baseline = "r00"),
#   df_temp_cons_c0 %>% filter(Type == "simulated") %>% mutate(Baseline = "c0"),
#   obs_temp_cons, "Temperature", "Algic", c(20,65)
# )


# # Save final dataset results for temp and NODF
# write.csv(obs_nodf_cons, "obs_nodf_cons_1000.csv", row.names = FALSE) 
# write.csv(obs_temp_cons, "obs_temp_cons_1000.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 8. CONSONANTS: Display the Plots ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Print combined plots for NODF and Temperature

ggsave('nodf_res_1000_cons_new_pval.png',
       plot=plot_nodf_cons,
       width=14,
       height=10,
       scale=0.9,
       bg='white')
# ggsave('temp_res_cons_1000.png',
#        plot=plot_temp_cons,
#        width=14,
#        height=10,
#        scale=0.9,
#        bg='white')

# 
# plot_temp_cons
plot_nodf_cons

