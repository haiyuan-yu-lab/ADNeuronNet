# High-Quality Violin Plot for Permutation Test Results
# Endophenotypes: Tau, Lysosome, Synapse, Mitochondrial
# Publication-ready figure with Nature-style aesthetics

setwd("./res/GSEA-Permutation/Pathobiological/All/")

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Define file names for the four endophenotypes of interest
file_configs <- list(
  list(
    file = "Tau_Fang-plus-Tracy-plus-GO_0.005.node_list",
    label = "Tau"
  ),
  list(
    file = "Lysosome_GO-plus-KEGG_0.005.node_list",
    label = "Lysosome"
  ),
  list(
    file = "Synaptic_GO-plus-KEGG_0.005.node_list",
    label = "Synapse"
  ),
  list(
    file = "Mitochondrial_Broad-plus-GO-plus-KEGG_0.005.node_list",
    label = "Mitochondrial"
  )
)

# Function to read and process each file
read_permutation_data <- function(config) {
  file_path <- paste0(config$file, "_permutation_results.csv")
  
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  
  data <- read.csv(file_path, comment.char = "#", stringsAsFactors = FALSE)
  
  # Extract observed odds ratio
  observed_or <- data %>%
    filter(query == "observed") %>%
    pull(Odds.Ratio)
  
  # Extract permuted odds ratios
  permuted_or <- data %>%
    filter(grepl("^perm_", query)) %>%
    pull(Odds.Ratio)
  
  # Calculate p-value
  p_value <- sum(permuted_or >= observed_or) / length(permuted_or)
  
  # Create data frame
  df <- data.frame(
    endophenotype = config$label,
    odds_ratio = permuted_or,
    type = "permuted"
  )
  
  # Add observed value
  observed_df <- data.frame(
    endophenotype = config$label,
    odds_ratio = observed_or,
    type = "observed"
  )
  
  list(
    permuted = df,
    observed = observed_df,
    p_value = p_value,
    label = config$label
  )
}

# Read all data
all_data <- lapply(file_configs, read_permutation_data)
all_data <- all_data[!sapply(all_data, is.null)]

if (length(all_data) == 0) {
  stop("No data files found. Please ensure the permutation result CSV files are in the working directory.")
}

# Combine permuted data
permuted_data <- bind_rows(lapply(all_data, function(x) x$permuted))

# Combine observed data
observed_data <- bind_rows(lapply(all_data, function(x) x$observed))

# Create p-value labels
p_labels <- data.frame(
  endophenotype = sapply(all_data, function(x) x$label),
  p_value = sapply(all_data, function(x) x$p_value)
)

# Set factor levels for consistent ordering
level_order <- c("Tau", "Lysosome", "Synapse", "Mitochondrial")
permuted_data$endophenotype <- factor(permuted_data$endophenotype, levels = level_order)
observed_data$endophenotype <- factor(observed_data$endophenotype, levels = level_order)
p_labels$endophenotype <- factor(p_labels$endophenotype, levels = level_order)

# Define color palette (Nature-style)
colors <- c(
  "Tau" = "#3C5488",
  "Lysosome" = "#00A087",
  "Synapse" = "#E64B35",
  "Mitochondrial" = "#F39B7F"
)

# Nature-style theme with legend text matching axis label size
theme_nature <- function() {
  theme_classic(base_size = 9, base_family = "sans") +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      axis.ticks.length = unit(0.12, "cm"),
      axis.title = element_text(size = 9, color = "black", face = "bold"),
      axis.text = element_text(size = 8, color = "black"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, face = "bold"),
      plot.margin = margin(10, 10, 10, 10, "pt"),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5)
    )
}

# Calculate y-axis limits and positions for annotations
y_max <- max(permuted_data$odds_ratio, observed_data$odds_ratio) * 1.18

# Create the violin plot
p <- ggplot(permuted_data, aes(x = endophenotype, y = odds_ratio, fill = endophenotype)) +
  # Violin with custom styling
  geom_violin(
    trim = FALSE,
    alpha = 0.7,
    color = "gray50",
    linewidth = 0.3,
    scale = "width",
    width = 0.75
  ) +
  # Add boxplot inside violin for summary statistics
  geom_boxplot(
    width = 0.12,
    fill = "white",
    color = "black",
    outlier.shape = NA,
    alpha = 0.9,
    linewidth = 0.4
  ) +
  # Add observed values as RED HORIZONTAL LINES (using geom_segment)
  geom_segment(
    data = observed_data,
    aes(
      x = as.numeric(endophenotype) - 0.15,
      xend = as.numeric(endophenotype) + 0.15,
      y = odds_ratio,
      yend = odds_ratio
    ),
    color = "#DC0000",
    linewidth = 1.2,
    inherit.aes = FALSE
  ) +
  # Add horizontal line at OR = 1 (no effect)
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "gray50",
    linewidth = 0.5
  ) +
  # Apply colors
  scale_fill_manual(values = colors) +
  # Labels
  labs(
    x = "Pathology",
    y = "Odds Ratio"
  ) +
  # Apply theme
  theme_nature() +
  # Adjust y-axis
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.1)),
    limits = c(NA, y_max)
  )

# Add p-value annotations above each violin
for (i in 1:nrow(p_labels)) {
  endo <- p_labels$endophenotype[i]
  pval <- p_labels$p_value[i]
  obs_or <- observed_data$odds_ratio[observed_data$endophenotype == endo]
  
  # Format p-value
  if (pval < 0.001) {
    p_text <- "P < 0.001"
  } else {
    p_text <- sprintf("P = %.3f", pval)
  }
  
  # Get max y for this group
  group_max <- max(permuted_data$odds_ratio[permuted_data$endophenotype == endo], obs_or)
  
  p <- p + annotate(
    "text",
    x = as.numeric(endo),
    y = group_max + (y_max - min(permuted_data$odds_ratio)) * 0.06,
    label = p_text,
    size = 2.5,
    fontface = "italic"
  )
}

# Add legend annotation for observed - TOP RIGHT position with RED LINE
# Legend text size = 9pt (same as axis label size)
p <- p + 
  # Red line in legend
  annotate(
    "segment",
    x = 3.55,
    xend = 3.85,
    y = y_max * 0.98,
    yend = y_max * 0.98,
    color = "#DC0000",
    linewidth = 1.2
  ) +
  # Legend text
  annotate(
    "text",
    x = 3.9,
    y = y_max * 0.98,
    label = "Observed",
    size = 9 / .pt,  # Convert pt to ggplot size (same as axis label)
    hjust = 0,
    color = "#DC0000",
    fontface = "bold"
  )

# Save outputs with updated dimensions (width=4, height=3)
ggsave("endophenotype_violin_plot.pdf", p, width = 4, height = 3, dpi = 300)
ggsave("endophenotype_violin_plot.png", p, width = 4, height = 3, dpi = 300)

cat("Done! Violin plot saved as:\n")
cat("  - endophenotype_violin_plot.pdf\n")
cat("  - endophenotype_violin_plot.png\n")