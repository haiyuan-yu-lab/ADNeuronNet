library(ggplot2)
library(dplyr)

# --- Configuration ---
PATH <- "/Users/houy2/Documents/projects/AD/AD_PPI/haiyuan_ppi_v6_Jan2026"
enrichr_dir <- "pathway_enrichr"
UpDate <- "Jan2026"

# --- Load Data ---
input_file <- file.path(PATH, enrichr_dir, 
                        paste0("ADNeuronNet_interactors_GOSt_results_filtered.csv"))
df <- read.csv(input_file, header = TRUE)
head(df)

# --- Prepare Data ---
# Use term_name directly
df$Pathway_Clean <- df$term_name

# Remove duplicates (keep the one with smallest adjusted p-value)
df <- df %>%
  group_by(Pathway_Clean) %>%
  arrange(adjusted_p_value) %>%
  slice(1) %>%
  ungroup()

# Calculate -log10 transformed p-value (already provided in the data)
df$neg_log10_pval <- df$negative_log10_of_adjusted_p_value

# Select top pathways
df_plot <- df %>% 
  filter(adjusted_p_value < 0.05) %>%
  arrange(adjusted_p_value) %>% 
  head(13)

# Calculate enrichment ratio (intersection_size / term_size)
df_plot$enrichment_ratio <- df_plot$intersection_size / df_plot$term_size

# --- Create Bubble Plot ---
p <- ggplot(df_plot, aes(x = neg_log10_pval, y = reorder(Pathway_Clean, neg_log10_pval))) +
  geom_point(aes(size = intersection_size, fill = enrichment_ratio), alpha = 0.7, shape = 21, color = "#0288D1", stroke = 0.5) +
  scale_size_continuous(name = "Overlap Genes Count", range = c(4, 12)) +
  scale_fill_gradient(low = "#80CBC4", high = "#0277BD", name = "Enrichment Ratio") +
  labs(title = "Pathway Enrichment Analysis",
       x = "-log10(Adjusted P-value)",
       y = NULL) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 9),
        legend.position = "right")

# --- Save Plot ---
output_file <- file.path(PATH, enrichr_dir, 
                         paste0("pathway_enrichment_HPA_", UpDate, ".pdf"))
ggsave(output_file, p, width = 6.5, height = 3.5)

print(paste("Plot saved to:", output_file))
print(paste("Total significant pathways (FDR < 0.05):", nrow(df_plot)))
