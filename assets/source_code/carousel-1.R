library(tidyverse)
library(RColorBrewer)

files <- list.files("../odp/step2-figures/ribbin_backup/", pattern = "rbh", full.names = T)

ortholog <- read_tsv(files[1]) |> 
  select(2,3)

for (file in files) {
  ortholog <- ortholog |> 
    full_join(read_tsv(file) |> select(2,3))
}

synteny <- read_tsv(files[str_detect(files, "Bowfin") & str_detect(files, "SpottedGar")] )
color <- synteny |>  count(Bowfin_scaf, SpottedGar_scaf) |>  arrange(desc(n))
color <- color |> 
  mutate(color = rep(brewer.pal(12, "Paired"), length.out = nrow(color)))
synteny <- select(synteny, -color) |> 
  left_join(select(color, -n))

ortholog <- ortholog |> 
  left_join(select(synteny, 2,3, color)) |> 
  mutate(color = case_when(
    is.na(color) ~ "#D3D3D3",
    T ~ color 
  ))

for (file in files) {
  rbh <- read_tsv(file) |> 
    select(-color) 
  
  cols <- colnames(ortholog)[colnames(ortholog) %in% colnames(rbh)]
  
  rbh <- rbh |> 
    left_join(
      select(ortholog, all_of(cols), color) |> na.omit() |> distinct()
    ) 
  
  print(file)
  print(count(rbh, color) |> arrange(desc(n)))
  
  write_tsv(rbh, (file |> str_replace("ribbin_backup", "synteny_nocolor")))
}
