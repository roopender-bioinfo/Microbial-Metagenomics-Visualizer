# Microbial-Metagenomics-Visualizer
This repository has the simple and easy script for the analysis of OTU and TAX files from dada2 or anyother pipeline. It was developed at EERL laboratory central university of rajasthan
# Load necessary libraries in your Rstudio
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(readxl)
library(fantaxtic)

# Load your OTU, TAX, and sample files (modify paths as necessary)
otu_table <- read.csv("path/to/otu_table.csv", row.names = 1)
tax_table <- read.csv("path/to/tax_table.csv", row.names = 1)
sample_data <- read.csv("path/to/sample_data.csv", row.names = 1)

# Convert to phyloseq objects
otu_ps <- otu_table(otu_table, taxa_are_rows = TRUE)
tax_ps <- tax_table(as.matrix(tax_table))
sample_ps <- sample_data(sample_data)

# Create phyloseq object
ps <- phyloseq(otu_ps, tax_ps, sample_ps)

# Alpha Diversity
alpha_div <- estimate_richness(ps, measures = c("Shannon", "Simpson"))
alpha_div_plot <- plot_richness(ps, x = "SampleType", measures = c("Shannon", "Simpson")) +
  geom_boxplot(aes(fill = SampleType)) +
  theme_minimal()

# Beta Diversity
ord <- ordinate(ps, method = "PCoA", distance = "bray")
beta_div_plot <- plot_ordination(ps, ord, color = "SampleType") +
  geom_point(size = 5) +
  theme_minimal()

# Function to sort and plot top taxa
plot_top_taxa <- function(ps, rank, top_n) {
  top_taxa <- names(sort(taxa_sums(ps), decreasing = TRUE)[1:top_n])
  ps_top <- prune_taxa(top_taxa, ps)
  bar_plot <- plot_bar(ps_top, fill = rank) +
    theme_minimal() +
    labs(title = paste("Top", top_n, rank))
  return(bar_plot)
}

# Plot top 10, 20, 30 taxa at different ranks
top_10_genus <- plot_top_taxa(ps, "Genus", 10)
top_20_family <- plot_top_taxa(ps, "Family", 20)
top_30_order <- plot_top_taxa(ps, "Order", 30)

# Nested bar plots using fantaxtic
nested_bar_plot <- fantaxtic_bar(ps, ranks = c("Phylum", "Class", "Order", "Family", "Genus"))

# Save plots
ggsave("alpha_diversity_plot.png", plot = alpha_div_plot)
ggsave("beta_diversity_plot.png", plot = beta_div_plot)
ggsave("top_10_genus_plot.png", plot = top_10_genus)
ggsave("top_20_family_plot.png", plot = top_20_family)
ggsave("top_30_order_plot.png", plot = top_30_order)
ggsave("nested_bar_plot.png", plot = nested_bar_plot)

# Print plots to the console
print(alpha_div_plot)
print(beta_div_plot)
print(top_10_genus)
print(top_20_family)
print(top_30_order)
print(nested_bar_plot)
