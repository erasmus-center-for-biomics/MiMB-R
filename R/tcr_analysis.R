# section 1

library(tidyverse)

if(!file.exists("data/assembled_data.rds"))
  stop("data/assembled_data.rds is not found, please run R/environment.R first")

if(!dir.exists("output"))
  dir.create("output")

denv <- readRDS('data/assembled_data.rds')


# section 2

filtered_tcr <- left_join(
  denv$tcr %>%
    group_by(cellid, locus) %>%
    mutate(
      rnk = rank(-reads, ties.method='first'),
      fraction = reads/sum(reads)) %>%
    ungroup() %>%
    filter(reads >= 10) %>%
    filter(rnk <= 2) %>%
    filter(productive) %>%
    filter(locus %in% c("TRA", "TRB")),
  denv$welllist %>% select(cellid, sample), 
  by = 'cellid')

write_tsv(filtered_tcr, "output/filtered_tcr.txt")


# section3

tcr_to_cells <- filtered_tcr %>%
  filter(rnk == 1) %>%
  group_by(sample, locus, v_call, d_call, j_call, cdr3_aa) %>%
  summarise(cells = n_distinct(cellid)) %>% ungroup() %>%
  group_by(sample, locus) %>%
  mutate(rnk = rank(-cells, ties.method = 'first')) %>% ungroup() 

write_tsv(tcr_to_cells, "output/tcr_to_cells.txt")

# section 4
tcrplt <- tcr_to_cells %>% 
  ggplot(aes(x=rnk, y=cells)) + 
  geom_point() + 
  facet_grid(sample ~ locus) +
  theme_light()
ggsave(tcrplt, filename="output/tcr_to_cells.png", width=6, height=8, dpi=300)