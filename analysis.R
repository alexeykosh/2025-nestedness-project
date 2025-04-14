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
library(cowplot)

theme_set(theme_bw())

# Setting main parameters
ALPHA_ <- 0.005
N_ITER_ <- 1000

#-------------------------------------------------------------------------------
# 1. Useful functions 
#-------------------------------------------------------------------------------

## Function to sort the matrix by columns and rows sum
sort_matrix <- function(M, iterations = 100) {
  for (i in seq_len(iterations)) {
    if (ncol(M) > 1) {
      M <- M[, order(colSums(M), decreasing = TRUE)]
    }
    if (nrow(M) > 1) {
      M <- M[order(rowSums(M), decreasing = TRUE), ]
    }
  }
  return(M)
}

## Function to display matrices
plot_matrix <- function(mat) {
  df <- melt(mat)
  colnames(df) <- c("Language", "Phoneme", "Value")
  # Keep the order of languages
  df$Language <- factor(df$Language, levels = rev(rownames(mat)))  
  ggplot(df, aes(x = Phoneme, y = Language, fill = factor(Value))) +
    geom_tile(color = "white") +
    scale_fill_manual(values = c("0" = "cornsilk", "1" = "black"), 
                      name = "Presence") +
    labs(x = "Phonemes", y = "Languages") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = 'None')
}

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
    ###############
    ###############
    # REMOVE FOR PROPER TESTING (after pre-reg)
    # df_one_hot_matrix <- matrix(rbinom(nrow_ * ncol_, 1, fill),
    #                             nrow = nrow_, ncol = ncol_,
    #                             dimnames = list(paste0("agent_",
    #                                                    LETTERS[1:nrow_]),
    #                                             paste0("item_", 1:ncol_)))
    ###############
    ###############
    print(fam_N)
    results <- oecosimu(df_one_hot_matrix, 
                        nestfun = function_type, 
                        method = shuffling_type,
                        parallel = -1,
                        nsimul = n_iter,
                        alternative = alt)
    sim <- results$oecosimu$simulated
    if (identical(function_type, nestednodf)) {
      new_rows <- data.frame(
        Family = rep(fam_N, length(sim)),
        Measure = rep('NODF', length(sim)),
        Type = rep('simulated', length(sim)),
        Value = as.numeric(sim),
        p_value = as.numeric(results$oecosimu$pval[2]),
        n_langs = nrow(df_one_hot_matrix)
      )
      new_rows_real <- data.frame(
        Family = fam_N,
        Measure = 'NODF',
        Type = 'real',
        Value = as.numeric(results$statistic$statistic[2]),
        p_value = as.numeric(results$oecosimu$pval[2]),
        n_langs = nrow(df_one_hot_matrix)
      )
    } else {
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
  
  # Create labels with colored bullets using Unicode and color codes
  axis_labels_df <- df_obs_final %>%
    distinct(Family, sig_cat) %>%
    mutate(
      bullet = case_when(
        sig_cat == "both" ~ "\u25CF",   # Black circle
        sig_cat == "one" ~ "\u25CF",     # Black circle
        sig_cat == "none" ~ "\u25CF"     # Black circle
      ),
      color_code = bullet_colors[sig_cat],
      label = paste0(Family, " ", bullet)
    )
  
  # Create named vector for labels
  labels_vector <- setNames(axis_labels_df$label, axis_labels_df$Family)
  
  ## Plot for baseline r00 (no legend)
  p_r00 <- ggplot(df_obs_final, aes(x = Value_r00, y = Family)) +
    geom_segment(
      data = sim_summary %>% filter(Baseline == "r00"),
      aes(x = lower, xend = upper, y = Family, yend = Family),
      color = "#A9A9A9", linewidth = 0.8
    ) +
    geom_point(
      aes(fill = shape_type_r00, shape = shape_type_r00),
      size = 3.5, color = "black"
    ) +
    scale_fill_manual(
      name = "Significance type",
      values = c("triangle" = "#64B5F6", "square" = "#FF6961", "circle" = "grey"),
      labels = c("Nested", "Antinested", "Not significant")
    ) +
    scale_shape_manual(
      name = "Significance type",
      values = c("triangle" = 24, "square" = 22, "circle" = 21),
      labels = c("Nested", "Antinested", "Not significant")
    ) +
    scale_y_discrete(labels = labels_vector) +
    scale_x_continuous(limits = c(0, 100)) +
    labs(x = x_label, y = "", title = "r00") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 10, color = axis_labels_df$color_code),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    )
  
  ## Plot for baseline c0 (will supply the legend)
  p_c0 <- ggplot(df_obs_final, aes(x = Value_c0, y = Family)) +
    geom_segment(
      data = sim_summary %>% filter(Baseline == "c0"),
      aes(x = lower, xend = upper, y = Family, yend = Family),
      color = "#A9A9A9", linewidth = 0.8
    ) +
    geom_point(
      aes(fill = shape_type_c0, shape = shape_type_c0),
      size = 3.5, color = "black"
    ) +
    scale_fill_manual(
      name = "Significance type",
      values = c("triangle" = "#64B5F6", "square" = "#FF6961", "circle" = "grey"),
      labels = c("Nested", "Antinested", "Not significant")
    ) +
    scale_shape_manual(
      name = "Significance type",
      values = c("triangle" = 24, "square" = 22, "circle" = 21),
      labels = c("Nested", "Antinested", "Not significant")
    ) +
    scale_x_continuous(limits = c(0, 100)) +
    scale_y_discrete(labels = NULL) +
    labs(x = x_label, title = "c0") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 9),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    )
  
  # Extract legend from p_c0 (the second plot)
  legend <- get_legend(
    p_c0 + 
      theme(legend.position = "bottom") +
      guides(fill = guide_legend(title = "Significance type"),
             shape = guide_legend(title = "Significance type"))
  )
  
  # Combine plots (without legend)
  plots_combined <- plot_grid(
    p_r00, p_c0,
    ncol = 2,
    align = "h",
    rel_widths = c(1.2, 1)
  )
  
  # Combine plots and legend
  final_plot <- plot_grid(
    plots_combined,
    legend,
    ncol = 1,
    rel_heights = c(1, 0.1)
  )
  
  return(final_plot)
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
  #########################################################################
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


#-------------------------------------------------------------------------------
# 2. Loading and processing data
#-------------------------------------------------------------------------------

# ## Read Phoible data
# df_values <- read.csv("data/cldf-datasets-phoible-f36deac/cldf/values.csv", 
#                       sep = ",", 
#                       header = TRUE)
# df_languages <- read.csv("data/cldf-datasets-phoible-f36deac/cldf/languages.csv", 
#                          sep = ",", 
#                          header = TRUE)
# ## Generate list of families 
# f_n <- df_languages %>% 
#   group_by(Family_Name) %>% 
#   summarise(n_l = n()) %>% 
#   filter(n_l > 3)  %>% 
#   filter(!Family_Name %in% c("", "Bookkeeping")) %>%
#   pull(Family_Name)
# 
# df_languages %>% 
#   group_by(Family_Name) %>% 
#   summarise(n_l = n()) %>% 
#   filter(n_l > 3)  %>%
#   filter(!Family_Name %in% c("", "Bookkeeping")) %>%
#   write.csv2(., file='data/summary.csv')

#-------------------------------------------------------------------------------
# 3. Run Nestedness Tests 
#-------------------------------------------------------------------------------

## Run NODF tests
# # For r00
# df_nodf_r00 <- nested_test(f_n, n_iter = N_ITER_, shuffling_type = "r00",
#                            function_type = nestednodf)
# # For c0
# df_nodf_c0  <- nested_test(f_n, n_iter = N_ITER_, shuffling_type = "c0",
#                            function_type = nestednodf)
# 
# ## Run Temperature tests
# # For r00
# df_temp_r00 <- nested_test(f_n, n_iter = N_ITER_, shuffling_type = "r00",
#                            function_type = nestedtemp)
# # For c0
# df_temp_c0  <- nested_test(f_n, n_iter = N_ITER_, shuffling_type = "c0",
#                            function_type = nestedtemp)

df_temp_c0 <- read.csv("df_temp_c0.csv")
df_temp_r00 <- read.csv("df_temp_r00.csv")
df_nodf_c0 <- read.csv("df_nodf_c0.csv")
df_nodf_r00 <- read.csv("df_nodf_r00.csv")

# # Save simulated datasets results for temp and NODF
# write.csv(df_temp_c0, "df_temp_c0.csv", row.names = FALSE)
# write.csv(df_temp_r00, "df_temp_r00.csv", row.names = FALSE)
# write.csv(df_nodf_c0, "df_nodf_c0.csv", row.names = FALSE)
# write.csv(df_nodf_r00, "df_nodf_r00.csv", row.names = FALSE)
# 
# # Save final dataset results for temp and NODF
# write.csv(obs_nodf, "obs_nodf.csv", row.names = FALSE)
# write.csv(obs_temp, "obs_temp.csv", row.names = FALSE)


# 

#-------------------------------------------------------------------------------
# 4. Dataset Manipulation and Plotting
#-------------------------------------------------------------------------------

## For NODF
# Process simulated and observed NODF data
df_nodf_r00_r <- df_nodf_r00 %>% mutate(significant = p_value <= ALPHA_)
df_nodf_c0_r  <- df_nodf_c0 %>% mutate(significant = p_value <= ALPHA_)
sim_nodf <- process_simulated(df_nodf_r00_r, df_nodf_c0_r)
obs_nodf <- process_real(df_nodf_r00_r, df_nodf_c0_r, "NODF", ALPHA_)
# Create combined NODF plot 
plot_nodf <- plot_combined(obs_nodf, sim_nodf$sim_summary, "NODF")
# Create NODF distribution plot
# dist_nodf <- plot_distribution(
#   df_nodf_r00 %>% filter(Type == "simulated") %>% mutate(Baseline = "r00"),
#   df_nodf_c0 %>% filter(Type == "simulated") %>% mutate(Baseline = "c0"),
#   obs_nodf, "NODF", "Algic", c(20, 65)
# )

## For Temperature
# Process simulated and observed Temperature data
df_temp_r00_r <- df_temp_r00 %>% mutate(significant = p_value <= ALPHA_)
df_temp_c0_r  <- df_temp_c0 %>% mutate(significant = p_value <= ALPHA_)
sim_temp <- process_simulated(df_temp_r00_r, df_temp_c0_r)
obs_temp <- process_real(df_temp_r00_r, df_temp_c0_r, "Temperature", ALPHA_)
# Create combined Temperature plot
plot_temp <- plot_combined(obs_temp, sim_temp$sim_summary, "Temperature")
# Create Temperature distribution plot
# dist_temp <- plot_distribution(
#   df_temp_r00 %>% filter(Type == "simulated") %>% mutate(Baseline = "r00"),
#   df_temp_c0 %>% filter(Type == "simulated") %>% mutate(Baseline = "c0"),
#   obs_temp, "Temperature", "Algic", c(20,65)
# )


#-------------------------------------------------------------------------------
# 5. Display the Plots
#-------------------------------------------------------------------------------

## Print combined plots for NODF and Temperature

ggsave('nodf_res_005.png', 
       plot=plot_nodf,
       width=14, 
       height=10,
       scale=0.9)
ggsave('temp_res_005.png', 
       plot=plot_temp,
       width=14,
       height=10,
       scale=0.9)



