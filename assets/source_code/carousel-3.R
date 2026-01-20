library(tidyverse)
library(readxl)
library(cowplot)
library(scales)
library(uwot)
library(Rtsne)
library(gridExtra)
library(grid)
library(pheatmap)

TGD_expression <- readRDS("../../datasets2/4.Expression_quantify/tpm.RDS")
SS4R_expression <- readRDS("../../datasets/5.Normalization/tpm.RDS")

species_list <- read_excel("../../table3.1.xlsx") |>
  select(-c(1:3)) |>
  mutate(latin = str_replace(`Scientific name`, " ", "_") |> tolower()) |>
  mutate(short = paste0(
    str_sub(`Scientific name`, 1, 1),
    str_extract(`Scientific name`, "(?<= ).{3}")
  ))

all_expression <- bind_rows(
  bind_rows(TGD_expression, .id = "spp") |>
    rownames_to_column("geneID") |>
    left_join(select(species_list, Species, latin),
      by = join_by(spp == latin)
    ) |>
    select(-spp),
  bind_rows(SS4R_expression[-c(1, 2)], .id = "spp") |>
    rownames_to_column("geneID") |>
    left_join(select(species_list, Species, short),
      by = join_by(spp == short)
    ) |>
    select(-spp)
)

# write_tsv(all_expression, "../../data/all_expression.tsv")


species_list <- species_list |>
  left_join(count(all_expression, Species, name = "Expressed")) |>
  drop_na() |>
  mutate(Expressed_rate = round(Expressed / `Longest isoforms`, 4)) |>
  mutate(Legend = paste0(Species, "\n", scales::comma(Expressed), " | ", scales::percent(Expressed_rate), "\n"))

p1 <- all_expression |>
  pivot_longer(bones:testis, names_to = "tissue", values_to = "tpm") |>
  mutate(
    tpm = log2(tpm + 1),
    Species = factor(Species, levels = species_list$Species)
  ) |>
  ggplot(aes(x = tpm, fill = Species, colour = Species)) +
  geom_density(alpha = 0.01, size = 0.2) +
  xlim(0.2, 7.5) +
  facet_wrap(~tissue, nrow = 5, scales = "free_y", strip.position = "bottom") +
  theme_cowplot() +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = "Log transformed transcript per million (log(2TPM+1))") +
  scale_fill_discrete(labels = species_list$Legend) +
  scale_color_discrete(labels = species_list$Legend) +
  theme(
    legend.spacing.y = unit(2, "cm"),
    legend.key.height = unit(1, "cm")
  )

ggsave("../pdf/Figure4.1.pdf", p1, width = 8, height = 10)

tissue_color <- c(
  "brain" = "#00a07b",
  "heart" = "#ca0000",
  "liver" = "#349a00",
  "testis" = "#ff6800",
  "bones" = "#2fcdfb",
  "kidney" = "#c99a00",
  "ovary" = "#be4c98",
  "gills" = "#6a329f",
  "intestine" = "#3297c7",
  "muscle" = "#ffcc00"
)

p2 <- all_expression |>
  pivot_longer(bones:testis, names_to = "tissue", values_to = "tpm") |>
  mutate(
    tpm = log2(tpm + 1),
    Species = factor(Species, levels = species_list$Species)
  ) |>
  ggplot(aes(x = tpm, fill = tissue, colour = tissue)) +
  geom_density(alpha = 0.01, size = 0.2) +
  xlim(0.2, 7.5) +
  facet_wrap(~Species, nrow = 5, scales = "free_y", strip.position = "bottom") +
  theme_cowplot() +
  scale_color_manual(values = tissue_color) +
  scale_fill_manual(values = tissue_color) +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.position = "top") +
  labs(x = "Log transformed transcript per million (log2TPM+1)")

ggsave("../appendix/appendix.fig_1.pdf", p2, width = 8, height = 10)


# UMAP ---------------------------------------------------------------------

species_shape <- c(
  "Spotted gar" = 16,
  "Bowfin" = 17,
  "Mexican tetra" = 15,
  "Atlantic cod" = 3,
  "European eel" = 4,
  "Zebrafish" = 8,
  "Medaka" = 7,
  "Northern pike" = 10,
  "Eastern mudminnow" = 12,
  "Grayling" = 2,
  "European whitefish" = 1,
  "American whitefish" = 5,
  "Brown trout" = 6,
  "Rainbow trout" = 13,
  "Brook trout" = 11
)

all_expression <- read_tsv("../../data/all_expression.tsv")
N2HOG <- read_tsv("../../datasets2/3.Ohnolog_pipeline/N2HOG.tsv")

# find the complete singleton
singleton_N2 <- N2HOG |>
  count(N2, spp) |>
  pivot_wider(names_from = spp, values_from = n, values_fill = 0) |>
  # select(-c(lepisosteus_oculatus, amia_calva )) |>
  filter(if_all(2:10, ~ .x == 1)) |>
  pull(N2)
# generate complete expression matrix
umap_input <- N2HOG |>
  filter(N2 %in% singleton_N2) |>
  # filter(!(spp %in% c("lepisosteus_oculatus", "amia_calva" ))) |>
  select(-category) |>
  left_join(all_expression) |>
  pivot_longer(bones:testis, names_to = "tissue") |>
  unite("spp_tissue", spp, tissue, sep = ".", remove = T) |>
  select(-geneID, -Species) |>
  pivot_wider(names_from = spp_tissue, values_from = value) |>
  column_to_rownames("N2") |>
  t()
# and log transform it
umap_input <- log2(umap_input + 1)

# apply umap
umap_result <- scale(umap_input) |>
  umap(n_neighbors = 15, min_dist = 0.05, n_components = 2)

umap_df <- as.data.frame(umap_result) |>
  rownames_to_column("spp_tissue") |>
  separate(spp_tissue, into = c("spp", "Tissue"), sep = "\\.", remove = T) |>
  left_join(select(species_list, latin, Species), by = join_by(spp == latin)) |>
  mutate(Species = factor(Species, levels = species_list$Species))

# UMAP Analysis of Tissue and Species Gene Expression Patterns

p3 <- ggplot(umap_df, aes(x = V1, y = V2, color = Tissue, shape = Species)) +
  geom_point(size = 3) +
  geom_text(aes(label = Tissue), size = 4, vjust = -1) +
  stat_ellipse(
    data = subset(umap_df, Species %in% c("Spotted gar", "Bowfin")),
    aes(group = Species, color = Species),
    level = 0.85, size = 1, linetype = "dashed"
  ) +
  geom_point(
    data = subset(umap_df, Species %in% c("Spotted gar", "Bowfin")),
    shape = 1, size = 6, color = "black", stroke = 1.5, alpha = 0.5
  ) +
  theme_cowplot() +
  scale_x_continuous(position = "top") +
  scale_color_manual(values = tissue_color) +
  scale_shape_manual(values = species_shape) +
  labs(
    title = "Telest WGD analysis",
    subtitle = "A",
    x = "UMAP1",
    y = "UMAP2"
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    # legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0)
  ) +
  guides(
    color = "none",
    shape = guide_legend(title = "", nrow = 2)
  ) +
  geom_text(
    data = subset(umap_df, Species == "Spotted gar"),
    aes(x = min(V1) * 0.9, y = max(V2) * 1.08, label = "Spotted gar"),
    inherit.aes = FALSE, color = "black", size = 4
  ) +
  geom_text(
    data = subset(umap_df, Species == "Bowfin"),
    aes(x = -0.2, y = max(V2) * 1.05, label = "Bowfin"),
    inherit.aes = FALSE, color = "black", size = 4, hjust = 0
  )


# tSNE on highly variable genes --------------------------------------------------------------------

variable_genes <- apply(umap_input, 2, var)
top_genes <- order(variable_genes, decreasing = TRUE)[1:1000]
tSNE_input <- umap_input[, top_genes]

tsne_result <- Rtsne(tSNE_input, dims = 2, perplexity = 10, theta = 0.5)
tsne_df <- as.data.frame(tsne_result$Y) |>
  setNames(c("tSNE1", "tSNE2")) |>
  cbind(select(umap_df, spp, Tissue)) |>
  left_join(select(species_list, latin, Species), by = join_by(spp == latin)) |>
  mutate(Species = factor(Species, levels = species_list$Species))


p4 <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = Tissue, shape = Species, group = Tissue)) +
  geom_point(size = 3) +
  geom_text(aes(label = Tissue), size = 4, vjust = -1) +
  geom_point(
    data = subset(tsne_df, Species %in% c("Spotted gar", "Bowfin")),
    shape = 1, size = 6, color = "black", stroke = 1.5, alpha = 0.5
  ) +
  theme_cowplot() +
  stat_ellipse(level = 0.95, aes(fill = Tissue), geom = "polygon", alpha = 0.1, color = NA) +
  scale_color_manual(values = tissue_color) +
  scale_fill_manual(values = tissue_color) +
  scale_shape_manual(values = species_shape) +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "grey95", color = NA)
  ) +
  labs(
    subtitle = "C",
    x = "t-SNE1",
    y = "t-SNE2"
  )


# Salmonid ----------------------------------------------------------------

N6HOG <- read_tsv("../../datasets/3.Ohnolog_pipeline/N6HOG.tsv")

singleton_N6 <- N6HOG |>
  count(N6, spp) |>
  pivot_wider(names_from = spp, values_from = n, values_fill = 0) |>
  # select(-c(esox_lucius , umbra_pygmae  )) |>
  filter(if_all(2:9, ~ .x == 1)) |>
  pull(N6)
# generate complete expression matrix
umap_input2 <- N6HOG |>
  filter(N6 %in% singleton_N6) |>
  # filter(!(spp %in% c("esox_lucius", "umbra_pygmae" ))) |>
  select(-category) |>
  left_join(all_expression) |>
  pivot_longer(bones:testis, names_to = "tissue") |>
  unite("spp_tissue", spp, tissue, sep = ".", remove = T) |>
  select(-geneID, -Species) |>
  pivot_wider(names_from = spp_tissue, values_from = value) |>
  column_to_rownames("N6") |>
  t()
# and log transform it
umap_input2 <- log2(umap_input2 + 1)

umap_result2 <- scale(umap_input2) |>
  umap(n_neighbors = 15, min_dist = 0.05, n_components = 2)

umap_df2 <- as.data.frame(umap_result2) |>
  rownames_to_column("spp_tissue") |>
  separate(spp_tissue, into = c("spp", "Tissue"), sep = "\\.", remove = T) |>
  left_join(select(species_list, latin, Species), by = join_by(spp == latin)) |>
  mutate(Species = factor(Species, levels = species_list$Species))

# MAP Analysis of Tissue and Species Gene Expression Patterns

p5 <-
  ggplot(umap_df2, aes(x = V1, y = V2, color = Tissue, shape = Species)) +
  geom_point(size = 3) +
  geom_text(aes(label = Tissue), size = 4, vjust = -1) +
  theme_cowplot() +
  scale_x_continuous(position = "top") +
  scale_y_continuous(position = "right") +
  scale_color_manual(values = tissue_color) +
  scale_shape_manual(values = species_shape) +
  labs(
    title = "Salmonid WGD analysis",
    subtitle = "B",
    x = "UMAP1",
    y = "UMAP2"
  ) +
  theme(
    legend.position = "bottom",
    legend.justification = c(1, 0),
    legend.box.just = "right",
    legend.box = "horizontal",
    # legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0)
  ) +
  guides(
    color = "none",
    shape = guide_legend(title = "", nrow = 2)
  ) +
  stat_ellipse(
    data = subset(umap_df2, Species %in% c("Northern pike", "Eastern mudminnow")),
    aes(group = Species, color = Species),
    level = 0.25, size = 1, linetype = "dashed"
  ) +
  geom_point(
    data = subset(umap_df2, Species %in% c("Northern pike", "Eastern mudminnow")),
    shape = 1, size = 6, color = "black", stroke = 1.5, alpha = 0.5
  ) +
  geom_text(
    data = subset(umap_df2, Species == "Northern pike"),
    aes(x = min(V1) * 0.75, y = mean(V2) * 2.25, label = "Northern pike"),
    inherit.aes = FALSE, color = "black", size = 4
  ) +
  geom_text(
    data = subset(umap_df2, Species == "Eastern mudminnow"),
    aes(x = mean(V1), y = max(V2) * 0.5, label = "Eastern\nmudminnow"),
    inherit.aes = FALSE, color = "black", size = 4, hjust = 0
  )

# tSNE on highly variable genes --------------------------------------------------------------------

variable_genes2 <- apply(umap_input2, 2, var)
top_genes2 <- order(variable_genes2, decreasing = TRUE)[1:1000]
tSNE_input2 <- umap_input2[, top_genes2]

tsne_result2 <- Rtsne(tSNE_input2, dims = 2, perplexity = 10, theta = 0.5)
tsne_df2 <- as.data.frame(tsne_result2$Y) |>
  setNames(c("tSNE1", "tSNE2")) |>
  cbind(select(umap_df2, spp, Tissue)) |>
  left_join(select(species_list, latin, Species), by = join_by(spp == latin)) |>
  mutate(Species = factor(Species, levels = species_list$Species))


p6 <-
  ggplot(tsne_df2, aes(x = tSNE1, y = tSNE2, color = Tissue, shape = Species, group = Tissue)) +
  geom_point(size = 3) +
  geom_text(aes(label = Tissue), size = 4, vjust = -1) +
  geom_point(
    data = subset(tsne_df2, Species %in% c("Northern pike", "Eastern mudminnow")),
    shape = 1, size = 6, color = "black", stroke = 1.5, alpha = 0.5
  ) +
  theme_cowplot() +
  stat_ellipse(level = 0.95, aes(fill = Tissue), geom = "polygon", alpha = 0.1, color = NA) +
  scale_y_continuous(position = "right") +
  scale_color_manual(values = tissue_color) +
  scale_fill_manual(values = tissue_color) +
  scale_shape_manual(values = species_shape) +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "grey95", color = NA)
  ) +
  labs(
    subtitle = "D",
    x = "t-SNE1",
    y = "t-SNE2"
  )


# final -------------------------------------------------------------------


final <- grid.arrange(p3, p5, p4, p6, ncol = 2, nrow = 2)
ggsave("../pdf/Figure4.2.pdf", final, width = 16, height = 16)
