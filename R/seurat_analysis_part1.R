
# section 1
library(tidyverse)
library(Seurat)


if(!file.exists("data/assembled_data.rds"))
  stop("data/assembled_data.rds is not found, please run R/environment.R first")

if(!dir.exists("output"))
  dir.create("output")

denv <- readRDS("data/assembled_data.rds")


# section 2
counts <- denv$expr %>%
  filter(gene_biotype == 'protein_coding') %>%
  pivot_wider(
    id_cols = gene_id, 
    names_from = 'cellid', 
    values_from = reads,
    values_fill = list('reads' = 0)) %>%
  (function(d) {
    m <- as.matrix(d[,2:ncol(d)])
    rownames(m) <- d$gene_id
    m
  })  


# section 3
scobj <- CreateSeuratObject(
    counts, project = "label", 
    min.cells = 1, min.features = 500)
scobj <- subset(scobj, cells = Cells(scobj))


# section 4 normalization and scaling
scobj <- NormalizeData(scobj, normalization.method = 'LogNormalize', scale.factor = 10000)
scobj <- ScaleData(scobj, features = rownames(scobj))
scobj <- FindVariableFeatures(scobj, selection.method='vst', nfeatures=2000)


# section 5 
scobj <- RunPCA(scobj, features=VariableFeatures(scobj))
scobj <- JackStraw(scobj, num.replicate=100)
scobj <- ScoreJackStraw(scobj, dims=1:20)


# section 6
ptbl <- tibble(
  component = paste0('PC', 1:length(scobj[['pca']]@stdev)),
  stdev = scobj[['pca']]@stdev) %>%
  mutate(
    component = parse_factor(component, levels = component),
    variance = stdev ^ 2,
    `variance %` = variance / sum(variance) * 100)

write_tsv(ptbl, "output/pca_stats.txt")
pplt <- ptbl %>%
  ggplot(aes(x=component, y=`variance %`)) +
  geom_bar(stat = 'identity') +
  theme_light() +
  theme(axis.text.x = element_text(angle=90))
ggsave(pplt, filename="output/pca_variance_explained.png", width=9, height=5, dpi=300)


jtbl <- as_tibble(scobj[['pca']]@jackstraw$overall.p.values) %>%
  mutate(
    component = parse_factor(
      str_c('PC', PC, sep = ''), levels=levels(ptbl$component)),
    `-log10(p-value)` = -log10(Score))
write_tsv(ptbl, "output/jackstraw_stats.txt")
jplt <- jtbl %>%
  ggplot(aes(x=component, y=`-log10(p-value)`)) +
  geom_bar(stat = 'identity') +
  theme_light() +
  theme(axis.text.x = element_text(angle=90))
ggsave(jplt, filename="output/jackstraw_score.png", width=9, height=5, dpi=300)


# save the seurat object
saveRDS(scobj,"output/seurat_obj.rds")