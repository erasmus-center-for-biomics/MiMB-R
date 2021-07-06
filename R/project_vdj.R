
library(tidyverse)
library(Seurat)


# check variables that need to be predefined
if(!"locus" %in% ls())
  stop("please specify the locus, TRA or TRB") 

if(!"coordinate_type" %in% ls())
  stop("please specify a coordinate_type, either tsne or umap")

if(coordinate_type == "tsne") {
  coordinate_file <- "output/tsne_coordinates.txt"
} else if(coordinate_type == "umap") {
  coordinate_file <- "output/umap_coordinates.txt"
}

# check files
if(!file.exists(coordinate_file))
  stop(sprintf("%s is not found, please run R/seurat_analysis_part2.R first", coordinate_file))

if(!file.exists("output/filtered_tcr.txt"))
  stop("output/filtered_tcr.txt is not found, please run R/tcr_analysis.R first")


if(!dir.exists("output"))
  dir.create("output")

# load the single cell and coordinate data
ftcr <- read_tsv("output/filtered_tcr.txt")
coordinates <- read_tsv(coordinate_file)


d <- left_join(
  coordinates,
  ftcr %>%
    filter(locus == !!locus) %>%
    filter(rnk == 1) %>%
    mutate(composition = sprintf("%s|%s|%s", v_call, d_call, j_call)),
  by = 'cellid')


tcrplt <- d %>%
  ggplot(aes(x=x, y=y)) +
  geom_point(size=0.5) +
  geom_point(aes(colour=composition), size=2, data=d%>% filter(!is.na(composition))) +
  theme_light() +
  theme(
    legend.position = 'None',
    aspect.ratio = 1,
    panel.grid = element_blank())

ggsave(
  tcrplt, 
  filename=file.path("output", sprintf("%s.%s.png", locus, coordinate_type)), 
  width=6, height=6, dpi=300)
  