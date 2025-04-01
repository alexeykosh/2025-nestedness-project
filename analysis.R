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

theme_set(theme_bw())

# Setting main parameters
ALPHA_ <- 0.005
N_ITER_ <- 100


# 1. Useful functions ---------------------------------------------------------

## Sort the matrix by columns and rows sum
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
      # axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      legend.position = 'None') 
}
## Function for nestendess testing
nested_test <- function(family_list, 
                        function_type = nestednodf, 
                        shuffling_type = 'r00',
                        save_ = TRUE,
                        n_iter = 2,
                        alt="two.sided"){
  # Create the empty dataframe
  df_results <- data.frame(Family = character(), 
                           Measure = character(), 
                           Type = character(),
                           Value = numeric(),
                           p_value = numeric(),
                           stringsAsFactors = TRUE)
  
  for (fam_N in family_list){
    ## Step 1: Select languages
    lgs <- df_languages[df_languages$Family_Name == fam_N,]$ID
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
          ## random tie-break if equal)
          min_count <- min(.x$n_phonemes)
          candidates <- .x %>% filter(n_phonemes == min_count)
          return(candidates %>% slice_sample(n = 1))
        } else {
          ## For more than 2 invetories: calculate the degree of overlap
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
    ## Convert to the matrix format
    df_one_hot_matrix <- as.matrix(df_one_hot[, -1])  # Remove Language_ID 
    rownames(df_one_hot_matrix) <- df_one_hot$Language_ID  # Set row names
    nrow_ <- nrow(df_one_hot_matrix)
    ncol_ <- ncol(df_one_hot_matrix)
    fill <- sum(df_one_hot_matrix) / (nrow_ * ncol_)
    ###############
    ###############
    # REMOVE FOR PROPER TESTING (after pre-reg)
    df_one_hot_matrix <- matrix(rbinom(nrow_ * ncol_, 1, fill),
                                nrow = nrow_, ncol = ncol_,
                                dimnames = list(paste0("agent_",
                                                       LETTERS[1:nrow_]),
                                                paste0("item_", 1:ncol_)))
    ###############
    ###############
    print(fam_N)
    results <- oecosimu(df_one_hot_matrix, 
                        nestfun=function_type, 
                        method=shuffling_type,
                        parallel=-1,
                        nsimul=n_iter,
                        alternative=alt)
    sim <- results$oecosimu$simulated
    if (identical(function_type, nestednodf)){
      new_rows <- data.frame(
        Family = rep(fam_N, length(sim)),
        Measure = rep('NODF', length(sim)),
        Type = rep('simulated', length(sim)),
        Value = as.numeric(sim),  # Ensure it's numeric
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
    }
    else{
      # if nested-temp
      new_rows <- data.frame(
        Family = rep(fam_N, length(sim)),
        Measure = rep('Temperature', length(sim)),
        Type = rep('simulated', length(sim)),
        Value = as.numeric(sim),  # Ensure it's numeric
        p_value = as.numeric(results$oecosimu$pval[2]),
        n_langs = nrow(df_one_hot_matrix)
      )
      new_rows_real <- data.frame(
        Family = fam_N,
        Measure = 'Temperature', 
        Type = 'real',
        Value = as.numeric(results$statistic$statistic[2]),
        p_value = as.numeric(results$oecosimu$pval[2]),
        n_langs = nrow(df_one_hot_matrix)
      )
    }
    df_results <- rbind(df_results, new_rows) 
    df_results <- rbind(df_results, new_rows_real) 
    
  }
  return(df_results)
}

# 2. Loading and processing data -----------------------------------------------
## Read Phoible data
df_values <- read.csv("data/cldf-datasets-phoible-f36deac/cldf/values.csv", 
                      sep=",", 
                      header=TRUE)
df_languages <- read.csv("data/cldf-datasets-phoible-f36deac/cldf/languages.csv", 
                         sep=",", 
                         header=TRUE)
## Generate list of families
f_n <- df_languages %>% 
  group_by(Family_Name) %>% 
  summarise(n_l=n()) %>% 
  filter(n_l > 3)  %>% 
  filter(! Family_Name %in% c("", "Bookkeeping")) %>%
  pull(Family_Name)

# 3. Results -----------------------------------------------------------------

## Get results
df_results_nodf_c0 <- nested_test(f_n, 
                                  n_iter=N_ITER_,
                                  shuffling_type = 'c0',
                                  function_type=nestednodf)
## Add singnificance
df_results_nodf_c0 <- df_results_nodf_c0 %>%
  mutate(significant = p_value <= ALPHA_) # Set to 0.005 for real data
## Save the data as csv
write.csv(df_results_nodf_c0,
          "df_results_nodf_c0_two_sided_small_family.csv",
          row.names = FALSE)
# Similarly for r00
df_results_nodf_r00 <- nested_test(f_n, 
                                  n_iter=N_ITER_,
                                  shuffling_type = 'r00',
                                  function_type=nestednodf)
## Add singnificance
df_results_nodf_r00 <- df_results_nodf_r00 %>%
  mutate(significant = p_value <= ALPHA_) # Set to 0.005 for real data
## Save the data as csv
write.csv(df_results_nodf_r00,
          "df_results_nodf_r00_two_sided_small_family.csv",
          row.names = FALSE)

# 4. Dataset manipulation for plots --------------------------------------------
## Load simulated data for r00 and c0
df_sim_r00 <- read.csv("df_results_nodf_r00_two_sided_small_family.csv") %>%
  filter(Type == "simulated") %>%
  mutate(Baseline = "r00")
df_sim_c0 <- read.csv("df_results_nodf_c0_two_sided_small_family.csv") %>%
  filter(Type == "simulated") %>%
  mutate(Baseline = "c0")
## Combine simulated datasets into one
df_sim <- bind_rows(df_sim_r00, df_sim_c0)
## Compute quantiles for simulated data
sim_summary <- df_sim %>%
  group_by(Family, Baseline) %>%
  summarise(lower = quantile(Value, 0.025),
            upper = quantile(Value, 0.975),
            .groups = "drop")
## Load and prepare real data for r00
df_obs_r00 <- read.csv("df_results_nodf_r00_two_sided_small_family.csv") %>%
  group_by(Family) %>%
  # Calculate simulated mean per family
  mutate(Simulated_Mean = mean(Value)) %>%
  # Determine the significant side
  mutate(Significant_Side = case_when(
    Type != "real" ~ NA_character_,
    p_value > ALPHA_ ~ NA,
    Value > Simulated_Mean ~ "Nested",
    Value < Simulated_Mean ~ "Antinested",
    TRUE ~ NA
  )) %>%
  filter(Type == "real") %>%
  # Assign a shape based on the direction of significance (Nested/Antinested)
  mutate(shape_type_r00 = case_when(
    Significant_Side == "Nested" ~ "triangle",   # Triangle for Nested
    Significant_Side == "Antinested" ~ "square", # Square for Antinested
    TRUE ~ "circle"                              # Circle if non significant
  )) %>%
  # Select the necessary columns
  select(Family, Value_r00 = Value, Significant_Side_r00 = Significant_Side, 
         Simulated_Mean_r00 = Simulated_Mean, significant_r00 = significant, 
         shape_type_r00, n_langs)
## Load and prepare real data for c0
df_obs_c0 <- read.csv("df_results_nodf_c0_two_sided_small_family.csv") %>%
  group_by(Family) %>%
  # Calculate simulated mean per family
  mutate(Simulated_Mean = mean(Value)) %>%
  # Determine the significant side
  mutate(Significant_Side = case_when(
    Type != "real" ~ NA_character_,
    p_value >= ALPHA_ ~ NA,
    Value > Simulated_Mean ~ "Nested",
    Value < Simulated_Mean ~ "Antinested",
    TRUE ~ NA
  )) %>%
  filter(Type == "real") %>%
  # Assign a shape based on the direction of significance (Nested/Antinested)
  mutate(shape_type_c0 = case_when(
    Significant_Side == "Nested" ~ "triangle",   # Triangle for Nested
    Significant_Side == "Antinested" ~ "square", # Square for Antinested
    TRUE ~ "circle"                              # Circle if non significant
  )) %>%
  # Select the necessary columns
  select(Family, Value_c0 = Value, Significant_Side_c0 = Significant_Side, 
         Simulated_Mean_c0 = Simulated_Mean, significant_c0 = significant, 
         shape_type_c0)
## Merge real data and create significance category
df_obs_final <- inner_join(df_obs_r00, df_obs_c0, by = "Family") %>%
  # See if a family is significant for none, one or both baselines
  mutate(sig_cat = case_when(
    significant_r00 & significant_c0 ~ "both",
    xor(significant_r00, significant_c0) ~ "one",
    TRUE ~ "none"
  ))
## Bullet point (with HTML for colors) for number of significant baselines
bullet_colors <- c("both" = "#009E73", "one" = "#FFA500", "none" = "black")
axis_labels_df <- df_obs_final %>% 
  distinct(Family, sig_cat) %>%
  mutate(label = paste0("<span style='color:black'>", Family, "</span> ",
                        "<span style='color:", 
                        bullet_colors[sig_cat], "'>&#9679;</span>"))
labels_vector <- setNames(axis_labels_df$label, axis_labels_df$Family)

# 5. Combined plot with Patchwork ----------------------------------------------
## Extract simulated summaries for each baseline
sim_summary_r00 <- sim_summary %>% filter(Baseline == "r00")
sim_summary_c0 <- sim_summary %>% filter(Baseline == "c0")
## Plot for baseline r00
p_r00 <- ggplot(df_obs_final, aes(x = Value_r00, y = Family)) +
  geom_segment(data = sim_summary_r00,
               aes(x = lower, xend = upper, y = Family, yend = Family),
               # dotted lines for simulated range
               color = "#A9A9A9", linewidth = 0.8) + 
  
  geom_point(aes(fill = shape_type_r00, shape = shape_type_r00),
             size = 3.5, color = "black") + # observed values
  
  ## Set the different colors
  scale_fill_manual(name = "Significance type",
                    values = c("triangle" = "#64B5F6", 
                               "square" = "#FF6961", 
                               "circle" = "grey"),
                    labels = c("triangle" = "Nested", 
                               "square" = "Antinested", 
                               "circle" = "Not significant")) +
  ## Set the different shapes
  scale_shape_manual(name = "Significance type",
                     values = c("triangle" = 24, 
                                "square" = 22, 
                                "circle" = 21),
                     labels = c("triangle" = "Nested", 
                                "square" = "Antinested", 
                                "circle" = "Not significant")) +
  
  scale_y_discrete(labels = labels_vector) +
  labs(x = "NODF", y = "Family", title = "r00") +
  scale_x_continuous(limits = c(0, NA)) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_markdown(size = 10),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")
## Plot for baseline c0
p_c0 <- ggplot(df_obs_final, aes(x = Value_c0, y = Family)) +
  geom_segment(data = sim_summary_c0,
               aes(x = lower, xend = upper, y = Family, yend = Family),
               # solid line for c0
               color = "#A9A9A9",  linewidth = 0.8) +
  geom_point(aes(fill = shape_type_c0, shape = shape_type_c0),
             size = 3.5, color = "black") +
  scale_fill_manual(name = "Significance type",
                    values = c("triangle" = "#64B5F6", 
                               "square" = "#FF6961", 
                               "circle" = "grey"),
                    labels = c("triangle" = "Nested", 
                               "square" = "Antinested", 
                               "circle" = "Not significant")) +
  scale_shape_manual(name = "Significance type",
                     values = c("triangle" = 24, 
                                "square" = 22, 
                                "circle" = 21),
                     labels = c("triangle" = "Nested", 
                                "square" = "Antinested", 
                                "circle" = "Not significant")) +
  # No legend for the y axis
  scale_y_discrete(labels = NULL) +
  # Set the axis and legend
  labs(x = "NODF", title = "c0") +
  scale_x_continuous(limits = c(0, NA)) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 9),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")
## Combine with patchwork
combined_plot <- (p_r00 + p_c0) +
  plot_layout(ncol = 2, guides = "auto") +
  plot_annotation(theme = theme(plot.title = element_text(hjust = 0.5))) &
  theme(legend.position = "bottom")
print(combined_plot)


# 6. Distribution figure for the entire dataset -------------------------

# ## Choose a family to preview 
# selected_family <- "Algic"  # (normally it's the entire phoible dataset)
# ## Filter simulated data for the selected family
# df_sim_r00_family <- df_sim_r00 %>% filter(Family == selected_family)
# df_sim_c0_family  <- df_sim_c0  %>% filter(Family == selected_family)
# ## Compute mean and standard deviation for each baseline
# mean_r00 <- mean(df_sim_r00_family$Value, na.rm = TRUE)
# sd_r00   <- sd(df_sim_r00_family$Value, na.rm = TRUE)
# mean_c0  <- mean(df_sim_c0_family$Value, na.rm = TRUE)
# sd_c0    <- sd(df_sim_c0_family$Value, na.rm = TRUE)
# ## Extract real values for the selected family
# real_r00 <- df_obs_r00 %>% filter(Family == selected_family) %>% pull(Value_r00)
# real_c0  <- df_obs_c0  %>% filter(Family == selected_family) %>% pull(Value_c0)
# ## Define x-axis range and density sequence
# x_lim <- c(15, 45)
# x_seq <- seq(0, 100, length.out = 200)
# ## Create data frames for the density curves
# df_density_r00 <- data.frame(
#   x = x_seq,
#   y = dnorm(x_seq, mean = mean_r00, sd = sd_r00),
#   baseline = "r00"
# )
# df_density_c0 <- data.frame(
#   x = x_seq,
#   y = dnorm(x_seq, mean = mean_c0, sd = sd_c0),
#   baseline = "c0"
# )
# ## Create distribution plot
# p_gauss <- ggplot() +
#   geom_area(data = df_density_r00, aes(x = x, y = y, fill = baseline), 
#             alpha = 0.3) +
#   geom_area(data = df_density_c0, aes(x = x, y = y, fill = baseline), 
#             alpha = 0.3) +
#   geom_line(data = df_density_r00, aes(x = x, y = y, color = baseline), 
#             linewidth = 1) +
#   geom_line(data = df_density_c0, aes(x = x, y = y, color = baseline), 
#             linewidth = 1) +
#   scale_fill_manual(values = c("r00" = "#56B", "c0" = "orange"), 
#                     guide = FALSE) +
#   scale_color_manual(name = "Baseline",
#                      values = c("r00" = "#56B", "c0" = "orange"),
#                      labels = c("r00", "c0")) +
#   labs(x = "NODF", y = "Density",
#        title = paste("Gaussian Distribution for the", selected_family, 
#                      "family")) +
#   scale_x_continuous(limits = x_lim) +
#   theme_minimal() +
#   theme(plot.title = element_text(hjust = 0.5))
# ## Add vertical line for the real value
# p_gauss <- p_gauss + 
#   geom_vline(xintercept = real_r00[1], color = "red", linetype = "dashed", 
#              linewidth = 1)
# print(p_gauss)