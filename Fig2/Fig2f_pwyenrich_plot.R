library(ggplot2)
library(dplyr)

# --- Configuration ---
PATH <- "/Users/houy2/Documents/projects/AD/AD_PPI/haiyuan_ppi_v6_Jan2026"
enrichr_dir <- "pathway_enrichr"
UpDate <- "Jan2026"

# --- Load Data ---
input_file <- file.path(PATH, enrichr_dir, 
                        paste0("PPI_pathway_enrichment_all_interactors_", UpDate, ".tsv"))
df <- read.delim(input_file, sep = "\t", header = TRUE)

# --- Prepare Data ---
# Remove pathway IDs (e.g., "hsa03051:") from pathway names
df$Pathway_Clean <- gsub("^hsa\\d+:", "", df$Pathway)

# Remove duplicates (keep the one with smallest adjusted p-value)
df <- df %>%
  group_by(Pathway_Clean) %>%
  arrange(Adjusted_P_Value) %>%
  slice(1) %>%
  ungroup()

# Calculate -log10 transformed p-value
df$neg_log10_pval <- -log10(df$Adjusted_P_Value)

# Select top pathways, excluding the specified pathway
df_plot <- df %>% 
  filter(Adjusted_P_Value < 0.05) %>%
  filter(Pathway_Clean != "ATP-dependent chromatin remodeling") %>%  # Remove specific pathway
  arrange(Adjusted_P_Value) %>% 
  head(10)  # Changed from 20 to 10

# --- Create Bubble Plot ---
p <- ggplot(df_plot, aes(x = neg_log10_pval, y = reorder(Pathway_Clean, neg_log10_pval))) +
  geom_point(aes(size = Overlap_Genes_Count, fill = Enrichment_Ratio), alpha = 0.7, shape = 21, color = "#E53935", stroke = 0.5) +
  scale_size_continuous(name = "Overlap Genes Count", range = c(4, 12)) +
  scale_fill_gradient(low = "#fecc5c", high = "#bd0026", name = "Enrichment Ratio") +
  labs(title = "Pathway Enrichment Analysis",
       x = "-log10(Adjusted P-value)",
       y = NULL) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 9),
        legend.position = "right")

# --- Save Plot ---
output_file <- file.path(PATH, enrichr_dir, 
                         paste0("pathway_enrichment_bubble_", UpDate, ".pdf"))
ggsave(output_file, p, width = 6.5, height = 3.5)

print(paste("Plot saved to:", output_file))
print(paste("Total significant pathways (FDR < 0.05):", nrow(df_plot)))
