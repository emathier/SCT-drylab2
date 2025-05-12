library(Seurat)
library(ggplot2)
library(tidyverse)

# load seurat obj
rm(list=ls())
seurat <- readRDS("~/projects/SCT-drylab2/Data/raw/seurat_for_qc.rds")

# Extract metadata
df <- seurat@meta.data %>% as_tibble(rownames = "cell")

# Free memory
rm(seurat)

# number of cells
cat("Number of cells",nrow(df))

# Convert to long format for ggplot
df <- df %>% pivot_longer(cols= c(nCount_RNA, nFeature_RNA, nCount_ATAC, nFeature_ATAC, percent.mt, nucleosome_signal, nucleosome_percentile, TSS.enrichment, TSS.percentile))

# Violin plot for number of features rna


# Calculate cell counts for each region
below_1000_count <- df %>%
  dplyr::filter(name == "nFeature_RNA", value < 1000) %>%
  nrow()

between_1000_8500_count <- df %>%
  dplyr::filter(name == "nFeature_RNA", value >= 1000, value <= 8500) %>%
  nrow()

above_8500_count <- df %>%
  dplyr::filter(name == "nFeature_RNA", value > 8500) %>%
  nrow()

# Create the violin plot
nrf <- ggplot(df %>% dplyr::filter(name == "nFeature_RNA"), aes(x = name ,y = value)) +
  geom_violin(fill = "#b2f5d1") + # Use a light mint color for the fill
  theme_minimal() + # Use a minimal theme for a clean look
  labs(
    x = "", # Correct x-axis label
    y = "", # y-axis label
    title = "Number of RNA Features" # Title of the plot
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16), # Center and size the title
    axis.title.x = element_text(size = 14), # Style for x-axis title
    axis.title.y = element_text(size = 14), # Style for y-axis title
    axis.text.x = element_text(size = 12),  # Style for x-axis text
    axis.text.y = element_text(size = 12),  # Style for y-axis text
  ) +
  # Add horizontal lines
  geom_hline(yintercept = 1000, color = "black", linetype = "dashed") +
  geom_hline(yintercept = 8500, color = "black", linetype = "dashed") +
  # Add colored areas
  annotate(
    "rect",
    ymin = -Inf,
    ymax = 1000,
    xmin = -Inf,
    xmax = Inf,
    fill = "red",
    alpha = 0.2 # Adjust transparency as needed
  ) +
  annotate(
    "rect",
    ymin = 8500,
    ymax = Inf,
    xmin = -Inf,
    xmax = Inf,
    fill = "red",
    alpha = 0.2 # Adjust transparency as needed
  ) +
  # Add annotations for cell counts
  annotate(
    "text",
    x = 1.4, # Adjust horizontal position as needed
    y = 500, # Adjust vertical position as needed
    label = paste("Cells < 1000:", below_1000_count),
    color = "black",
    size = 4
  ) +
  annotate(
    "text",
    x = 1.4, # Adjust horizontal position as needed
    y = 5000, # Adjust vertical position as needed
    label = paste("Cells 1000-8500:", between_1000_8500_count),
    color = "black",
    size = 4
  ) +
  annotate(
    "text",
    x = 1.4, # Adjust horizontal position as needed
    y = 10000, # Adjust vertical position as needed
    label = paste("Cells > 8500:", above_8500_count),
    color = "black",
    size = 4
  )
  
ggsave(filename = "~/projects/SCT-drylab2/Deliverables/polishedQC/nFeaturesRNA.png", units = "mm" , width = 150, height = 150, dpi  = 400)


# Percent mitochondiral
# Assuming your data frame is named 'df'
# And it has columns 'name' and 'value'

# Calculate cell counts for each region
below_20_count <- df %>%
  dplyr::filter(name == "percent.mt", value < 20) %>%
  nrow()

above_20_count <- df %>%
  dplyr::filter(name == "percent.mt", value > 20) %>%
  nrow()

# Create the violin plot
nrf <- ggplot(df %>% dplyr::filter(name == "percent.mt"), aes(x = name ,y = value)) +
  geom_violin(fill = "#b2f5d1") + # Use a light mint color for the fill
  theme_minimal() + # Use a minimal theme for a clean look
  labs(
    x = "", # Correct x-axis label
    y = "", # y-axis label
    title = "Percent Mitochondrial Reads" # Title of the plot
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16), # Center and size the title
    axis.title.x = element_text(size = 14), # Style for x-axis title
    axis.title.y = element_text(size = 14), # Style for y-axis title
    axis.text.x = element_text(size = 12),  # Style for x-axis text
    axis.text.y = element_text(size = 12),  # Style for y-axis text
  ) +
  # Add horizontal line
  geom_hline(yintercept = 20, color = "black", linetype = "dashed") +
  # Add colored area
  annotate(
    "rect",
    ymin = 20,
    ymax = Inf,
    xmin = -Inf,
    xmax = Inf,
    fill = "red",
    alpha = 0.2 # Adjust transparency as needed
  ) +
  # Add annotations for cell counts
  annotate(
    "text",
    x = 1.4, # Adjust horizontal position as needed
    y = 10, # Adjust vertical position as needed
    label = paste("Cells < 20%:", below_20_count),
    color = "black",
    size = 4
  ) +
  annotate(
    "text",
    x = 1.4, # Adjust horizontal position as needed
    y = 30, # Adjust vertical position as needed
    label = paste("Cells > 20%:", above_20_count),
    color = "black",
    size = 4
  )

ggsave(filename = "~/projects/SCT-drylab2/Deliverables/polishedQC/precentmt.png", units = "mm" , width = 150, height = 150, dpi  = 400)




# Calculate cell counts for each region
below_1000_count_atac <- df %>%
  dplyr::filter(name == "nFeature_ATAC", value < 1000) %>%
  nrow()

between_1000_30000_count_atac <- df %>%
  dplyr::filter(name == "nFeature_ATAC", value >= 1000, value <= 30000) %>%
  nrow()

above_30000_count_atac <- df %>%
  dplyr::filter(name == "nFeature_ATAC", value > 30000) %>%
  nrow()

# Create the violin plot
atac_plot <- ggplot(df %>% dplyr::filter(name == "nFeature_ATAC"), aes(x = name ,y = value)) +
  geom_violin(fill = "#b2f5d1") + # Use a light mint color for the fill
  theme_minimal() + # Use a minimal theme for a clean look
  labs(
    x = "",
    y = "",
    title = "Number of Features ATAC"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
  ) +
  # Add horizontal lines
  geom_hline(yintercept = 1000, color = "black", linetype = "dashed") +
  geom_hline(yintercept = 30000, color = "black", linetype = "dashed") +
  # Add colored areas
  annotate(
    "rect",
    ymin = -Inf,
    ymax = 1000,
    xmin = -Inf,
    xmax = Inf,
    fill = "red",
    alpha = 0.2
  ) +
  annotate(
    "rect",
    ymin = 30000,
    ymax = Inf,
    xmin = -Inf,
    xmax = Inf,
    fill = "red",
    alpha = 0.2
  ) +
  # Add annotations for cell counts
  annotate(
    "text",
    x = 1.4,
    y = 500,
    label = paste("Cells < 1000:", below_1000_count_atac),
    color = "black",
    size = 4
  ) +
  annotate(
    "text",
    x = 1.4,
    y = 15000,
    label = paste("Cells 1000-30000:", between_1000_30000_count_atac),
    color = "black",
    size = 4
  ) +
  annotate(
    "text",
    x = 1.4,
    y = 35000,
    label = paste("Cells > 30000:", above_30000_count_atac),
    color = "black",
    size = 4
  )

print(atac_plot)

ggsave(filename = "~/projects/SCT-drylab2/Deliverables/polishedQC/nFeaturesATAC.png", units = "mm" , width = 150, height = 150, dpi  = 400)


# Violin plot for tss enrichment

# Assuming your data frame is named 'df'
# And it has columns 'name' and 'value'

# Calculate cell counts for each region
below_1_count_tss <- df %>%
  dplyr::filter(name == "TSS.enrichment", value < 1) %>%
  nrow()

above_1_count_tss <- df %>%
  dplyr::filter(name == "TSS.enrichment", value > 1) %>%
  nrow()

# Create the violin plot
tss_plot <- ggplot(df %>% dplyr::filter(name == "TSS.enrichment"), aes(x = name ,y = value)) +
  geom_violin(fill = "#b2f5d1") + # Use a light mint color for the fill
  theme_minimal() + # Use a minimal theme for a clean look
  labs(
    x = "",
    y = "",
    title = "TSS Enrichment"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
  ) +
  # Add horizontal line
  geom_hline(yintercept = 1, color = "black", linetype = "dashed") +
  # Add colored area *below* the line
  annotate(
    "rect",
    ymin = -Inf,  # Changed from 1 to -Inf
    ymax = 1,      # Changed from Inf to 1
    xmin = -Inf,
    xmax = Inf,
    fill = "red",
    alpha = 0.2
  ) +
  # Add annotations for cell counts
  annotate(
    "text",
    x = 1.4,
    y = 0.5,
    label = paste("Cells < 1:", below_1_count_tss),
    color = "black",
    size = 4
  ) +
  annotate(
    "text",
    x = 1.4,
    y = 2,
    label = paste("Cells > 1:", above_1_count_tss),
    color = "black",
    size = 4
  )

print(tss_plot)

ggsave(filename = "~/projects/SCT-drylab2/Deliverables/polishedQC/TSSEnrichment.png", units = "mm" , width = 150, height = 150, dpi  = 400)

library(ggplot2)
library(dplyr)

# Assuming your data frame is named 'df'
# And it has columns 'name' and 'value'

# Calculate cell counts for each region
below_1_count_ns <- df %>%
  dplyr::filter(name == "nucleosome_signal", value < 1) %>%
  nrow()

above_1_count_ns <- df %>%
  dplyr::filter(name == "nucleosome_signal", value > 1) %>%
  nrow()

# Create the violin plot
ns_plot <- ggplot(df %>% dplyr::filter(name == "nucleosome_signal"), aes(x = name ,y = value)) +
  geom_violin(fill = "#b2f5d1") + # Use a light mint color for the fill
  theme_minimal() + # Use a minimal theme for a clean look
  labs(
    x = "",
    y = "Nucleosome Signal",
    title = "Nucleosome Signal"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
  ) +
  # Add horizontal line
  geom_hline(yintercept = 1, color = "black", linetype = "dashed") +
  # Add colored area *above* the line
  annotate(
    "rect",
    ymin = 1,      # Changed from -Inf to 1
    ymax = Inf,    # Changed from 1 to Inf
    xmin = -Inf,
    xmax = Inf,
    fill = "red",
    alpha = 0.2
  ) +
  # Add annotations for cell counts
  annotate(
    "text",
    x = 1.4,
    y = 0.5,
    label = paste("Cells < 1:", below_1_count_ns),
    color = "black",
    size = 4
  ) +
  annotate(
    "text",
    x = 1.4,
    y = 2,
    label = paste("Cells > 1:", above_1_count_ns),
    color = "black",
    size = 4
  )

print(ns_plot)

ggsave(filename = "~/projects/SCT-drylab2/Deliverables/polishedQC/NucleosomeSignal.png", units = "mm" , width = 150, height = 150, dpi  = 400)
