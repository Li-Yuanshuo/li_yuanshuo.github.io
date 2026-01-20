library(tidyverse)
library(ape)
library(ggtree)
library(pheatmap)
library(scales)
library(gridExtra)

# load datasets -----------------------------------------------------------
species_list <- read_csv("../../datasets2/telest_species.csv")

species_list <- species_list |>
  mutate(scientific_name = scientific_name |>
    tolower() |>
    str_replace(" ", "_"))



# F1.A tree ---------------------------------------------------------------
tree <- read.tree("../../datasets2/3.Ohnolog_pipeline/wkdir/OrthoFinder/Results_Dec11/Gene_Duplication_Events/SpeciesTree_Gene_Duplications_0.5_Support_figure.txt")
plot(tree, show.node.label = TRUE)
identify(tree)

subtree <- tree |>
  drop.tip(c(1,2))
plot(subtree)

tips_order <- species_list |>
  arrange(
    str_remove(subtree$tip.label, ".{5}$") |>
      match(species_list$scientific_name |> str_replace(" ", "_") |> tolower())
  ) |> 
  pull(common_name)  
  
tips_order[9] <- "Eastern\nmudminnow"

node_label <- tree$node.label[-2] |>
  str_replace("_", ":") 
node_label[1] <- "N0:3140\n"
node_label[2] <- "N2:1454\n"
node_label[4] <- "N4:1450\n"

p1 <- ggtree(tree |> 
               drop.tip(c(1)))  +
   geom_tiplab(aes(label = c(c("Mammalia", tips_order), rep("NA", 9))), align = T, size = 4.5, hjust = 0) +
  
  geom_nodelab(aes(label = c(rep("NA", 10), node_label)), hjust = -0.05, size = 3.5) +
  expand_limits(x = max(tree$edge.length) * 3) +
  scale_y_reverse() +
  geom_cladelabel(node = 13, label = "Closest pre-WGD", align = T, barsize = 1, fontsize = 3, 
                  offset = 0.2, angle = 90, hjust = 0.5, offset.text = 0.02, geom = "label") +
    geom_cladelabel(node = 14, label = "post-WGD species", align = T, barsize = 1, fontsize = 3, 
                    offset = 0.2, angle = 90, hjust = 0.5, offset.text = 0.02, geom = "label") 

TGD_node <- ggplot_build(p1)$data[[1]] |> filter(node == 14)


p1 <- p1 +
  geom_label(aes(x = (TGD_node$x + TGD_node$xend) / 2, y = -TGD_node$y),
    label = "TGD",
    fill = "orange",
    color = "black",
    label.size = 0.6,
    label.padding = unit(0.3, "lines"),
    label.r = unit(0.15, "lines"),
    hjust = 0.5,
    size = 3
  ) +
  geom_treescale(
    x = 0.05, y = -8.5, width = 0.1, fontsize = 3, linesize = 0.5,
    label = "substitutions\nper site", offset = 0.05,
    offset.label = -0.1
  ) +
  theme_tree() +
  theme(text = element_text(lineheight = 0.3))


p1 <- p1 +
  ggtitle("A") 
ggsave("../subfigures/figure3.1_A.png", p1, width = 5, height = 6, dpi = 300)

# F1.B heatmap ------------------------------------------------------------

ohno_list <- read_tsv("../../datasets2/3.Ohnolog_pipeline/dataset_final.tsv")

expressed_genes <- list.files("../../datasets2/4.Expression_quantify/expressed_gene_list/",
                              full.names = T) |>
  map_dfr(~ read_tsv(.x, show_col_types = FALSE))

OGtbl <- read_tsv("../../datasets2/3.Ohnolog_pipeline/wkdir/OrthoFinder/OGtbl.tsv", show_col_types = FALSE) |> 
  filter(geneID %in% expressed_genes$gene) |> 
  filter(!is.na(N2)) 

# select(OGtbl, geneID, spp, N2) |> 
#   left_join(select(ohno_list, N2HOG, category), by = join_by(N2 == N2HOG)) |> 
#   write_tsv("../../datasets2/3.Ohnolog_pipeline/N2HOG.tsv")

OGs_statis <- OGtbl |> 
  count(N2, spp) |> 
  pivot_wider(names_from = spp, values_from = n) %>%
  mutate_all(~ replace(., is.na(.), 0)) |> 
  left_join(select(ohno_list, N2HOG, category), by = join_by(N2 == N2HOG)) |> 
  mutate(OG_type = case_when(
    category %in% c("post-spec/wgd", "pre-spec/wgd") ~ "ohnolog",
    category == "singleton" ~ "singleton",
    T ~ "other"
  ))
           

species_order <- species_list |>
  arrange(c(2, 1, 3, 4, 5, 8, 9 , 6, 7)) |>
  pull(2) |> str_replace(" ", "_") |> tolower()

heatmap_input <- OGs_statis |>
  column_to_rownames("N2") |>
  select(all_of(species_order)) %>%
  mutate(across(everything(), ~ ifelse(. >= 3, 3, .))) |>
  t()

annotation_col <- OGs_statis |>
  column_to_rownames("N2") |>
  select(OG_type)

annotation_colors <- list(
  OG_type = c(
    "other" = "white",
    "ohnolog" = "#e76325",
    #"ohnolog_partial" = "#ffc069",
    "singleton" = "#5f3e98"
    #"singleton_partial" = "#b2abd3"
  )
)

p2 <- pheatmap(heatmap_input,
  color = c("white", "#b3b4b4", "#4c4c4d", "#040707"),
  breaks = c(-0.1, 0.5, 1.5, 2.5, 3.5),
  cluster_rows = F,
  cluster_cols = T,
  scale = "none",
  show_rownames = F,
  show_colnames = F,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  legend = F,
  annotation_legend = F,
  gaps_row = 2
)

ggsave("../subfigures/figure3.1_B.png", p2, width = 6, height = 6, dpi = 300)

# F1.C gene count ---------------------------------------------------------

No.genes <- read_tsv("../../datasets2/3.Ohnolog_pipeline/wkdir/OrthoFinder/Results_Dec11/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv",
  n_max = 3
) |>
  rename(number = 1) |>
  select(number, species_order) |>
  slice(3) |>
  mutate(number = "not in orthogroup") |>
  bind_rows(OGs_statis |>
    group_by(OG_type) |>
    summarise(across(where(is.numeric), sum, na.rm = TRUE)) |>
    rename(number = OG_type) |>
    select(number, species_order))

p3 <- No.genes |>
  pivot_longer(cols = -number, names_to = "species", values_to = "gene_count") |>
  mutate(species = factor(species, levels = rev(species_order))) |>
  mutate(number = factor(number |> 
                           str_replace("not in orthogroup",
                                       "not in \northogroup"), 
                         levels = rev(c("singleton", "ohnolog", "other", "not in \northogroup")))) |> 
  
  ggplot(aes(x = species, y = gene_count, fill = number)) +
  geom_col(color = "black", width = 0.7) +
  coord_flip() +
  labs(x = NULL, y = "Numbers of genes", fill = NULL) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    legend.position = c(1, 1.15), 
    legend.justification = c(1, 1), 
    legend.box.margin = margin(0, 0, 0, 0),
    plot.margin = margin(50, 0, 0, 0)
  ) +
  scale_y_continuous(labels = label_number(scale = 0.001, suffix = "k")) +
  scale_fill_manual(values = c(
    "not in \northogroup" = "black",
    "other" = "white",
    "ohnolog" = "#e76325",
    #"ohnolog_partial" = "#ffc069",
    "singleton" = "#5f3e98"
    #"singleton_partial" = "#b2abd3"
  )) 

ggsave("../subfigures/figure3.1_C.png", p3, width = 4, height = 6, dpi = 300)
