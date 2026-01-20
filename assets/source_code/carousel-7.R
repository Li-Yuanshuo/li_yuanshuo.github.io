library(tidyverse)
library(readxl)
library(cowplot)
library(scales)
library(uwot)
library(Rtsne)
library(gridExtra)
library(grid)
library(pheatmap)

# heatmap -----------------------------------------------------------------

species_color <- c(
  "Spotted gar" = "#4A90E2", "Bowfin" = "#005BBB", 
  "Mexican tetra" = "#F54E42", "Atlantic cod" = "#E63B2E", "European eel" = "#D32727",
  "Zebrafish" = "#C41D1D", "Medaka" = "#B01212", "Northern pike" = "#9A0808",
  "Eastern mudminnow" = "#840101", # 
  "Grayling" = "#6ABF69", "European whitefish" = "#4DAE4C", "American whitefish" = "#32A032",
  "Brown trout" = "#198519", "Rainbow trout" = "#106B10", "Brook trout" = "#085508" # 
)

ann_colors <- list(
  Tissue = tissue_color,
  Species = species_color
)


# TGD
row_annotation <- data.frame(
  Tissue = gsub(".*\\.(.*)", "\\1", rownames(umap_input)),
  latin = gsub("(.*)\\..*", "\\1", rownames(umap_input))
) |>
  left_join(select(species_list, latin, Species)) |>
  select(-latin)

rownames(row_annotation) <- rownames(umap_input)

heatmap1 <- pheatmap(
  umap_input,
  scale = "row",
  useRaster = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  annotation_row = row_annotation,
  annotation_colors = ann_colors,
  show_rownames = F,
  show_colnames = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(50)
)
library(ggplotify)

ggsave("../subfigures/figure4.0_A.png", as.ggplot(heatmap1) +
  ggtitle("A   Telest WGD analysis"), width = 12, height = 8, dpi = 300)

# ss4R
row_annotation2 <- data.frame(
  Tissue = gsub(".*\\.(.*)", "\\1", rownames(umap_input2)),
  latin = gsub("(.*)\\..*", "\\1", rownames(umap_input2))
) |>
  left_join(select(species_list, latin, Species)) |>
  select(-latin)

rownames(row_annotation2) <- rownames(umap_input2)

heatmap2 <- pheatmap(
  umap_input2,
  scale = "row",
  useRaster = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  annotation_row = row_annotation2,
  annotation_colors = ann_colors,
  show_rownames = F,
  show_colnames = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  legend = FALSE,
  annotation_legend = FALSE
)

ggsave("../subfigures/figure4.0_B.png", as.ggplot(heatmap2) +
  ggtitle("B   Salmonid WGD analysis"), width = 9.5, height = 8, dpi = 300)
