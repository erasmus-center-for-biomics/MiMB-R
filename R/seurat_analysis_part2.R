
#
if(! "ndim" %in% ls()) 
  stop("ndim not set")

if(! "scobj" %in% ls()) 
  stop("Please run seurat_analysis_part1.R first")


# section 7
coordinates <- RunUMAP(scobj, dims = 1:ndim) %>%
  (function(x){
    a <- as.data.frame(x[['umap']]@cell.embeddings)
    colnames(a) <- c('x', 'y')
    left_join(
      as_tibble(a) %>% mutate(cellid = rownames(a)),
      denv$welllist, by = 'cellid')
  })

write_tsv(coordinates, "output/umap_coordinates.txt")

umapplt <- coordinates %>% 
  ggplot(aes(x=x, y=y, colour=sample)) + 
  geom_point() + 
  theme_light() +
  theme(
    panel.grid = element_blank(),
    aspect.ratio = 1)
ggsave(umapplt, filename="output/umap.png", width=6, height=5, dpi=300)


# section 7
coordinates <- RunTSNE(scobj, dims = 1:ndim) %>%
  (function(x){
    a <- as.data.frame(x[['tsne']]@cell.embeddings)
    colnames(a) <- c('x', 'y')
    left_join(
      as_tibble(a) %>% mutate(cellid = rownames(a)),
      denv$welllist, by = 'cellid')
  })

write_tsv(coordinates, "output/tsne_coordinates.txt")

tsneplt <- coordinates %>% 
  ggplot(aes(x=x, y=y, colour=sample)) + 
  geom_point() + 
  theme_light() +
  theme(
    panel.grid = element_blank(),
    aspect.ratio = 1)
ggsave(tsneplt, filename="output/tsne.png", width=6, height=5, dpi=300)
