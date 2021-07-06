# section 1

library(tidyverse) 

denv <- new.env() 

denv$welllist <- read_tsv("data/welllist.txt") %>% 
  select(sample=Sample, row=Row, col=Col, barcode=Barcode, source_well=`Source well`) %>% 
  mutate(cellid = sprintf("%s:%s:%d:%d", sample, barcode, row, col)) 


# section 2
report_files <- file.path("data/TCR/", list.files("data/TCR/", pattern=".report.csv")) 
report_files <- report_files[!str_detect(report_files, '_NA_NA')]
matched <- str_match(report_files, "_([0-9]+)_([0-9]+).report.csv") 
rows <- parse_integer(matched[,2])
columns <- parse_integer(matched[,3])


# section 3
denv$tcr <- left_join( 
  1:length(report_files) %>% 
    map_df(function(idx){ 
      tryCatch({
        aln <- read_delim(report_files[idx], ';', col_types = 'ccclcccd')
        if(nrow(aln) > 0){
          aln %>%
            group_by(v_call, d_call, j_call, productive, locus) %>%
            summarise(
              reads = sum(n),
              top = max(n) / reads,
              n_cdr3_aa = n_distinct(cdr3_aa),
              cdr3_aa = cdr3_aa[which.max(n)[1]]
            ) %>%
            ungroup() %>%
            arrange(desc(reads)) %>%
            mutate(row=rows[idx], col=columns[idx])
        } else {
          NULL
        }
      }, error = function(e){NULL})
    }),
  denv$welllist %>% select(row, col, cellid),
  by = c("row", "col")
) %>% filter(!is.na(cellid))

# section 4
header <- c('type', 'gene_biotype', 'gene_id', 'symbol', 'name', 'barcode', 'row', 'col', 'reads')
denv$expr <- left_join(
  file.path("data/Expression/", list.files("data/Expression/", pattern="exon.tsv.gz")) %>%
  map_df(function(path){
    read_tsv(path, col_names = header, col_types = 'cccccciii') %>%
      select(-row, -col)
  }),
  denv$welllist %>% select(barcode, cellid),
  by = c("barcode")
) %>% filter(!is.na(cellid))

# section 5
saveRDS(denv, "data/assembled_data.rds")

