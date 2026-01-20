library(tidyverse)
library(readxl)
library(cowplot)
library(ggpubr)

species_list <- read_excel("../../table3.1.xlsx") |>
  select(-c(1:3)) |>
  mutate(group = c(rep("Mammalia", 2), rep("Holostei", 2), rep("Teleost", 7), rep("Salmonid", 6)) |>
    factor(levels = c("Mammalia", "Holostei", "Teleost", "Salmonid"))) |> 
  mutate(Species = factor(Species, levels = c("Human",  "Mouse",
                                              "Bowfin", "Spotted gar",
                                              "European eel", "Zebrafish", "Mexican tetra", "Atlantic cod", "Medaka","Northern pike", "Eastern mudminnow",
                                              "Grayling", "European whitefish", "American whitefish", "Brook trout", "Rainbow trout", "Brown trout"
                                              )))

Expression <- read_tsv("../../data/all_expression.tsv")

calculate_tau <- function(values) {
  if (all(values == 0)) {
    return(NA)
  }
  max_value <- max(values)
  normalized_values <- (1 - (values / max_value))
  tau <- sum(normalized_values) / (length(values) - 1)
  return(tau)
}

Tau <- Expression %>%
  mutate(across(bones:testis, ~ log2(. + 1), .names = "log2_{.col}")) %>%
  rowwise() %>%
  mutate(tau = calculate_tau(c_across(starts_with("log2_")))) %>%
  ungroup() |>
  select(geneID, Species, tau)

write_tsv(Tau, "../../data/Tau.tsv")
Tau <- read_tsv("../../data/Tau.tsv")
#  p1 ---------------------------------------------------------------------


p1_data <- Tau |>
  left_join(select(species_list, Species, group)) |> 
  mutate(Species = factor(Species, levels = c("Human",  "Mouse",
                                              "Bowfin", "Spotted gar",
                                              "European eel", "Zebrafish", "Mexican tetra", "Atlantic cod", "Medaka","Northern pike", "Eastern mudminnow",
                                              "Grayling", "European whitefish", "American whitefish", "Brook trout", "Rainbow trout", "Brown trout"
  )))


p1 <- p1_data |>
  ggplot(aes(x = tau, colour = Species, fill = Species)) +
  geom_density(aes(size = ifelse(Species %in% c("Northern pike", "Eastern mudminnow"), 1.5, 0.5)), alpha = 0) +
  scale_size_continuous(range = c(0.5, 1.5), guide = "none") +
  geom_boxplot(aes(y = 0.3),
    width = 0.8, outlier.shape = NA, alpha = 0.3,
    position = position_dodge2(preserve = "single")
  ) +
  facet_wrap(~group) +
  theme_cowplot() +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.position = "none") +
  labs(
    x = expression("Tissue Specificity (" * tau * ")"),
    y = "density"
  ) +
  ggtitle("A")


label_positon <- data.frame(
  Species = species_list$Species[-c(1:2)],
  position = ggplot_build(p1)$data[[2]]$y
) |>
  mutate(group = c(rep("Holostei", 2), rep("Teleost", 7), rep("Salmonid", 6)) |>
    factor(levels = c("Mammalia", "Holostei", "Teleost", "Salmonid")))

p1 <- p1 + geom_text(
  data = label_positon,
  aes(x = 0.05, y = position, label = Species),
  inherit.aes = FALSE,
  size = 3.5, hjust = 0, vjust = -0.5
)


ggsave("../subfigures/figure4.4_A.png", p1, width = 14, height = 5, dpi = 300)


# p2 ortholog Tau between pre-WGD species----------------------------------------------------------------------

# TGD ---------------------------------------------------------------------

N2HOG <- read_tsv("../../datasets2/3.Ohnolog_pipeline/N2HOG.tsv")

N2_tau <- Tau |>
  inner_join(select(N2HOG, -spp)) |>
  na.omit() |>
  filter(category != "pre-spec/ssd") |>
  mutate(category = case_when(
    str_detect(category, "wgd") ~ "ohnolog",
    T ~ category
  ))

p2.1_data <- N2_tau |>
  filter(Species %in% c("Spotted gar", "Bowfin")) |>
  select(-geneID) |>
  group_by(N2, Species) |>
  filter(n() <= 1) |>
  group_by(N2) |>
  filter(n_distinct(Species) == 2) |>
  ungroup() |>
  pivot_wider(names_from = Species, values_from = tau, values_fill = NA) 
  

regression_results1 <- p2.1_data %>%
  group_by(category) %>%
  summarise(
    slope = coef(lm(`Spotted gar` ~ Bowfin))[2],
    r_squared = summary(lm(`Spotted gar` ~ Bowfin))$r.squared,
    .groups = "drop"
  )

p2.1_data <- p2.1_data |>
  left_join(regression_results1, by = "category")


p2.1 <- p2.1_data |>
  ggplot(aes(x = Bowfin, y = `Spotted gar`, colour = category)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, size = 1) +
  scale_color_manual(values = c("ohnolog" = "#37b4f9", "singleton" = "#e69f00")) +
  scale_fill_manual(values = c("ohnolog" = "#37b4f9", "singleton" = "#e69f00")) +
  labs(
    title = "B          TGD Analysis",
    x = expression("Bowfin " * tau),
    y = expression("Spotted gar " * tau)
  ) +
  stat_regline_equation(
    aes(label = ..eq.label.., colour = category),
    size = 4,
    show.legend = FALSE
  ) +
  stat_cor(
    aes(label = ..r.label.., colour = category),
    size = 4, hjust = -3,
    show.legend = FALSE
  ) +
  theme_cowplot() +
  theme(
    legend.position = "bottom",
    legend.key = element_blank(),
    legend.text = element_text(size = 16)
  ) +
  guides(
    colour = guide_legend(title = NULL, override.aes = list(size = 10, shape = 16))
  )

# SS4R

N6HOG <- read_tsv("../../datasets/3.Ohnolog_pipeline/N6HOG.tsv")

N6_tau <- Tau |>
  inner_join(select(N6HOG, -spp)) |>
  na.omit() |>
  filter(category != "pre-spec/ssd") |>
  mutate(category = case_when(
    str_detect(category, "wgd") ~ "ohnolog",
    T ~ category
  ))

p2.2_data <- N6_tau |>
  filter(Species %in% c("Northern pike", "Eastern mudminnow")) |>
  select(-geneID) |>
  group_by(N6, Species) |>
  filter(n() <= 1) |>
  group_by(N6) |>
  filter(n_distinct(Species) == 2) |>
  ungroup() |>
  pivot_wider(names_from = Species, values_from = tau, values_fill = NA)

regression_results2 <- p2.2_data %>%
  group_by(category) %>%
  summarise(
    slope = coef(lm(`Northern pike` ~ `Eastern mudminnow`))[2],
    r_squared = summary(lm(`Northern pike` ~ `Eastern mudminnow`))$r.squared,
    .groups = "drop"
  )

p2.2_data <- p2.2_data |>
  left_join(regression_results1, by = "category")


p2.2 <- p2.2_data |>
  ggplot(aes(x = `Northern pike`, y = `Eastern mudminnow`, colour = category)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, size = 1) +
  scale_color_manual(values = c("ohnolog" = "#37b4f9", "singleton" = "#e69f00")) +
  scale_fill_manual(values = c("ohnolog" = "#37b4f9", "singleton" = "#e69f00")) +
  labs(
    title = "C          SS4R Analysis",
    x = expression("Northern pike " * tau),
    y = expression("Eastern mudminnow " * tau)
  ) +
  stat_regline_equation(
    aes(label = ..eq.label.., colour = category),
    size = 4,
    show.legend = FALSE
  ) +
  stat_cor(
    aes(label = ..r.label.., colour = category),
    size = 4, hjust = -3,
    show.legend = FALSE
  ) +
  theme_cowplot() +
  theme(
    legend.position = "none",
    legend.key = element_blank(),
    legend.text = element_text(size = 10)
  ) +
  guides(
    colour = guide_legend(title = NULL, override.aes = list(size = 5, shape = 16))
  )

p2 <- plot_grid(p2.1, p2.2, ncol = 1, align = "v")
ggsave("../subfigures/figure4.4_BC.png", p2, width = 8, height = 16, dpi = 300)



# p3 ----------------------------------------------------------------------

#TGD
N2_concerved_tau <- p2.1_data |> 
  mutate(tau_abs = abs(Bowfin - `Spotted gar`),
         tau_avg = (Bowfin +`Spotted gar`)/2) |> 
  filter(tau_abs < 0.1) |> 
  select(N2, category, tau_avg)

p3.1_data <- N2_tau |> 
  filter(!(Species %in% c("Spotted gar", "Bowfin"))) |> 
  inner_join(N2_concerved_tau) |> 
  mutate(tau_change = tau - tau_avg) |> 
  mutate(Species = factor(Species, levels = c("Human",  "Mouse",
                                              "Bowfin", "Spotted gar",
                                              "European eel", "Zebrafish", "Mexican tetra", "Atlantic cod", "Medaka","Northern pike", "Eastern mudminnow",
                                              "Grayling", "European whitefish", "American whitefish", "Brook trout", "Rainbow trout", "Brown trout"
  )))


p3.1 <- p3.1_data |>
  ggplot(aes(x = category, y = tau_change, color = category)) +
  geom_boxplot(aes(fill = category), alpha = 0.5, outlier.shape = NA, color = "grey30") +
  geom_violin(aes(fill = category), size = 1, alpha = 0.05, trim = TRUE) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 0.8) +
  stat_compare_means(method = "wilcox.test",size = 2.5) +
  ylim(-0.5, 0.5) +
  facet_wrap(~Species, nrow = 1) +
  scale_color_manual(values = c("ohnolog" = "#37b4f9", "singleton" = "#e69f00")) +
  scale_fill_manual(values = c("ohnolog" = "#37b4f9", "singleton" = "#e69f00")) +
  theme_cowplot() +
  theme(
    legend.position = "none", 
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(),  
  ) +
  labs(y = expression("decrease"%<-%"             "*tau*" change after WGD               "%->%"increase")) +
  ggtitle("D")
  
#SS4R

N6_concerved_tau <- p2.2_data |> 
  mutate(tau_abs = abs(`Northern pike` - `Eastern mudminnow`),
         tau_avg = (`Northern pike` + `Eastern mudminnow`)/2) |> 
  filter(tau_abs < 0.1) |> 
  select(N6, category, tau_avg)

p3.2_data <- N6_tau |> 
  filter(!(Species %in% c("Northern pike", "Eastern mudminnow"))) |> 
  inner_join(N6_concerved_tau) |> 
  mutate(tau_change = tau - tau_avg) |> 
  mutate(Species = factor(Species, levels = c("Human",  "Mouse",
                                              "Bowfin", "Spotted gar",
                                              "European eel", "Zebrafish", "Mexican tetra", "Atlantic cod", "Medaka","Northern pike", "Eastern mudminnow",
                                              "Grayling", "European whitefish", "American whitefish", "Brook trout", "Rainbow trout", "Brown trout"
  )))


p3.2 <- p3.2_data |>
  ggplot(aes(x = category, y = tau_change, color = category)) +
  geom_boxplot(aes(fill = category), alpha = 0.5, outlier.shape = NA, color = "grey30") +
  geom_violin(aes(fill = category), size = 1, alpha = 0.05, trim = TRUE) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 0.8) +
  stat_compare_means(method = "wilcox.test", size = 3) +
  ylim(-0.5, 0.5) +
  facet_wrap(~Species, nrow = 1) +
  scale_color_manual(values = c("ohnolog" = "#37b4f9", "singleton" = "#e69f00")) +
  scale_fill_manual(values = c("ohnolog" = "#37b4f9", "singleton" = "#e69f00")) +
  theme_cowplot() +
  theme(
    legend.position = "none", 
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(),  
  ) +
  labs(y = expression("decrease"%<-%"             "*tau*" change after WGD               "%->%"increase")) +
  ggtitle("E")
 
p3 <- plot_grid(p3.1, p3.2, ncol = 1, align = "v")
ggsave("../subfigures/figure4.4_DE.png", p3, width = 16, height = 10, dpi = 300)


plot_grid(p1,
          plot_grid(p2, p3, ncol = 2, rel_widths = c(1, 2)), 
          ncol = 1, rel_heights = c(1, 2.5))

ggsave("../pdf/Figure4.4.pdf", width = 16, height = 16, dpi = 300)

