
library(tidyverse)
library(scales)
library(Seurat)
library(wesanderson)

# check variables that need to be predefined
if(!"gene" %in% ls())
  stop("please specify a gene")

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

if(!file.exists("output/seurat_obj.rds"))
  stop("output/seurat_obj.rds is not found, please run R/seurat_analysis_part1.R first")

if(!dir.exists("output"))
  dir.create("output")

# load the single cell and coordinate data
scobj <- readRDS("output/seurat_obj.rds")
coordinates <- read_tsv(coordinate_file)


# get the data for gene
mat <- as.matrix(GetAssayData(scobj, 'scale.data'))
tf <- rownames(mat) == gene

if(sum(tf) == 0)
  stop(sprintf("Gene %s not found", gene))


d <- left_join(
  coordinates,
  mat[tf,] %>%
    (function(v, g=gene){
      tibble(
        cellid=names(v),
        gene=g,
        scaled_expression=v,
        expression = rescale(v, c(0,1)))
    }),
  by = 'cellid')


geneplt <- d %>%
  ggplot(aes(x=x, y=y, colour=expression)) +
  geom_point() +
  scale_colour_gradientn(colours = wes_palette('Zissou1', 5, 'continuous')) +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    aspect.ratio = 1)

ggsave(
  geneplt, 
  filename=file.path("output", sprintf("%s.%s.png", gene, coordinate_type)), 
  width=6, height=6, dpi=300)
  