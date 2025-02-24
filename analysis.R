# Vegan Package Description
# P.51 : null model (r1, c0, curveball, etc.)
# P. 145 : nested functions
# P. 153 : oecosimu (p value + simulations)

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

theme_set(theme_bw())


# Functions ---------------------------------------------------------------

## 1. Sort the matrix by columns and rows sum
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

## 2. Function to display matrices
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


# Loading data ------------------------------------------------------------

## Read Phoible
df_values <- read.csv("data/cldf-datasets-phoible-f36deac/cldf/values.csv", 
                      sep=",", 
                      header=TRUE)
df_languages <- read.csv("data/cldf-datasets-phoible-f36deac/cldf/languages.csv", 
                         sep=",", 
                         header=TRUE)


# Analysis -----------------------------------------------------------

## Generate list of families
f_n <- df_languages %>% 
  group_by(Family_Name) %>% 
  summarise(n_l=n()) %>% 
  # filter(n_l > 3)  %>%
  # filter(n_l < 10)  %>%
  filter(n_l > 3)  %>% 
  filter(! Family_Name %in% c("", "Bookkeeping")) %>%
  pull(Family_Name)

## Function for nestendess testing
nested_test <- function(family_list, 
                        function_type = nestednodf, 
                        shuffling_type = 'r00',
                        save_ = TRUE,
                        n_iter = 2){
  
  df_results <- data.frame(Family = character(), 
                           Measure = character(), 
                           Type = character(),
                           Value = numeric(),
                           p_value = numeric(),
                           stringsAsFactors = TRUE)
  
  for (fam_N in family_list){
    
    lgs <- df_languages[df_languages$Family_Name == fam_N,]$ID
    # lgs <- df_languages$ID
    
    ## Step 2: Select random (Language_ID, Contribution_ID) pairs
    selected_combinations <- df_values %>%
      dplyr::filter(Language_ID %in% lgs) %>%
      dplyr::group_by(Language_ID) %>%
      dplyr::slice_sample(n = 1) %>%  
      dplyr::ungroup() %>%
      dplyr::select(Language_ID, Contribution_ID)
    
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
    
    ## Convert to matrix
    df_one_hot_matrix <- as.matrix(df_one_hot[, -1])  # Remove Language_ID 
    rownames(df_one_hot_matrix) <- df_one_hot$Language_ID  # Set row names
    
    nrow_ <- nrow(df_one_hot_matrix)
    ncol_ <- ncol(df_one_hot_matrix)
    
    fill <- sum(df_one_hot_matrix) / (nrow_ * ncol_)
    
    # REMOVE FOR PROPER TESTING (after pre-reg)
    df_one_hot_matrix <- matrix(rbinom(nrow_ * ncol_, 1, fill),
                                nrow = nrow_, ncol = ncol_,
                                dimnames = list(paste0("agent_", LETTERS[1:nrow_]),
                                                paste0("item_", 1:ncol_)))
    
    print(fam_N)
    # print(df_one_hot_matrix)
    
    if (identical(function_type, nestednodf)){
      # bigger NODF -- nestedness
      # smaller NODF -- no nestedness
      alt = 'greater'
    }
    else{
      # smaller temperature (0) -- greater nestedness, 
      # bigger temp (100) -- no nestendess
      alt = 'less'
    }
    
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
        p_value = as.numeric(results$oecosimu$pval[3])
      )
      
      new_rows_real <- data.frame(
        Family = fam_N,
        Measure = 'NODF',
        Type = 'real',
        Value = as.numeric(results$statistic$statistic[3]),
        p_value = as.numeric(results$oecosimu$pval[3])
      )
      
    }
    else{
      # if nested-temp
      new_rows <- data.frame(
        Family = rep(fam_N, length(sim)),
        Measure = rep('NODF', length(sim)),
        Type = rep('simulated', length(sim)),
        Value = as.numeric(sim),  # Ensure it's numeric
        p_value = as.numeric(results$oecosimu$pval[1])
      )
      
      new_rows_real <- data.frame(
        Family = fam_N,
        Measure = 'NODF',
        Type = 'real',
        Value = as.numeric(results$statistic$statistic[1]),
        p_value = as.numeric(results$oecosimu$pval[1])
      )
    }
    
    df_results <- rbind(df_results, new_rows) 
    df_results <- rbind(df_results, new_rows_real) 
    
  }
  
  return(df_results)
  
}

# Results -----------------------------------------------------------------

## Get results
df_results_nodf_r00 <- nested_test(f_n, 
                                   n_iter=100,
                                   shuffling_type = 'c0',
                                   function_type=nestedtemp)
## Add singnificance
df_results_nodf_r00 <- df_results_nodf_r00 %>%
  mutate(significant = p_value < 0.005)

## Create a named vector for Family label colors
family_colors <- df_results_nodf_r00 %>%
  distinct(Family, significant) %>%
  mutate(color = ifelse(significant, "red", "black")) %>%
  { setNames(.$color, .$Family) }  # Use setNames instead of deframe()

## Figure opt. 1: color of the y-labels based on significance. 
ggplot(df_results_nodf_r00 %>% filter(Type == 'simulated'), 
       aes(x = Value, y = Family)) +
  stat_pointinterval() +  
  xlim(0, 100) +
  geom_point(data = df_results_nodf_r00 %>% 
               filter(Type == 'real') %>% 
               distinct(Family, Value), 
             aes(x = Value, y = Family), 
             color = "blue", size = 3) +  
  theme_minimal() +
  labs(y = '',
       x = 'NODF') +
  scale_y_discrete(labels = function(f) {
    # Apply color formatting to y-axis labels
    sapply(f, function(label) {
      if (label %in% names(family_colors)) {
        paste0("<span style='color:", family_colors[label], 
               "'>", label, "</span>")
      } else {
        label
      }
    })
  }) +
  theme(axis.text.y = ggtext::element_markdown())  

## Figure opt. 2: color of the points based on p value 
# Define color palette
pal <- wes_palette("Zissou1", 100, type = "continuous")

# Create plot
ggplot(df_results_nodf_r00 %>% filter(Type == 'simulated'), 
       aes(x = Value, y = Family)) +
  stat_pointinterval() +  
  xlim(0, 100) +
  geom_point(data = df_results_nodf_r00 %>% 
               filter(Type == 'real') %>% 
               distinct(Family, Value, p_value), 
             aes(x = Value, y = Family, color = p_value), 
             size = 3) +  
  theme_minimal() +
  labs(y = '',
       x = 'Temperature') +
  theme(legend.position = "bottom") +
  # scale_color_gradient2(low = "red", midpoint = 0.05, 
  #                       mid = "yellow", high = "blue") +
  scale_color_gradientn(colors = pal) 


