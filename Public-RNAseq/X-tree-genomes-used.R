# Query entrez for public RNA-seq data

library(tidyverse)
library(ape)
library(ggtree)

library(xml2)
library(rentrez)

library(patchwork)
library(ggpubr)

################################################################################
# genomes used after filtering etc

dat <-
  '3-RNA-Browser/*txid*' |>
  Sys.glob() |>
  basename() |>
  as_tibble() |>
  separate(value, c('name', 'txid', 'bio', 'order'), sep = '_') |>
  mutate_at('txid', str_remove, 'txid')

################################################################################
# highlight in species tree?


ref.trees <-
  '~/OCyRS-pipeline/reference-trees/*.tree' |>
  Sys.glob() %>%
  set_names(. %>% basename %>% fs::path_ext_remove()) |>
  map(read.tree)


ref.trees |>
  map('tip.label') |>
  map(length) |>
  unlist()
# Chen_genomic Chen_multigene         Cornet          Moore 
# 487            487            364            100 


# compare tax in trees to those in rnaseq
x.tax <-
  ref.trees |>
  map('tip.label') |>
  map(str_remove, '\\..*$') |>
  map(unique)
  

x.tax$`RNA-seq` <-
  dat |>
  pull(txid)


p.venn <- venn::venn(x.tax, zcolor = 'style',
                     ilabels = 'counts',
                     ilcs = 1, sncs = 1.2,
                     ggplot =  TRUE)

dat |>
  filter(txid == setdiff(x.tax$`RNA-seq`, x.tax$Chen_genomic))
# name                                 txid  bio          order
# <chr>                                <chr> <chr>        <chr>
# 1 Prochlorococcus.marinus.str.MIT.9312 74546 SAMN02598321 NA   
# Not part in the chen genomic tree

################################################################################
# location of species with rnaseq for selection in species tree

ref <- ref.trees$Chen_genomic
ref$tip.label <- str_remove(ref$tip.label, '\\..*$')

# ncbi order for coloration to accept 
rec.orders <-
  '../OCyRS-pipeline/data/A_representatives/taxonomy.tsv' |>
  read_tsv() |>
  select(order) |>
  drop_na() |>
  unique() |>
  mutate_at('order', str_remove, '^[0-9]* ')


# query all ranks for taxid
rank.helper <- function(x) {
  # x <- ncbi.res[[1]]
  # x |>
  #   read_xml() |>
  #   xml_child() -> x
  tibble(
    tax = x |>
      xml_child() |>
      xml_text(),
    ranks = x|>
      xml_find_all(".//*[name()='ScientificName']") |>
      map(xml_text) |>
      unlist()
  )
}

# query all clade names
ncbi.res <-
  ref$tip.label %>%
  split(floor(seq_along(.) / 50))  |>
  map(~ entrez_fetch('taxonomy', id = .x, rettype = 'xml'))

ncbi.ranks <-
  ncbi.res |>
  map(read_xml) |>
  map(xml_children) |>
  invoke(.f = c) |>
  map(rank.helper) |>
  bind_rows()


# filter for order
ncbi.orders <-
  ncbi.ranks |>
  semi_join(rec.orders, c('ranks' = 'order')) |>
  dplyr::rename(order = ranks) |>
  # make structure for groupOTU
  select(tax, order) |>
  unique() |>
  # add entry without order info
  complete(tax = ref$tip.label) |>
  mutate_at('order', replace_na, 'unknown') |>
  group_by(order) |>
  do(i = list(.$tax)) |>
  with(set_names(i, order)) |>
  map(1)

# for adding rnaseq info
rna.extra <-
  '0-public-stat.tsv' |>
  read_tsv() |>
  filter(no.libs >= 5) |>
  mutate_at('tax', as.character) |>
  left_join(
    dat |>
      transmute(txid, with.rnaseq = 'QC passed'),
    c('tax' = 'txid')
  ) |>
  mutate_at('with.rnaseq', replace_na, 'QC failed')



p.tree <- ggtree(ref, layout = 'circular', branch.length = 'none') %<+% rna.extra


p.otu <-
  groupOTU(p.tree, ncbi.orders, 'Order') +
  aes(color = Order) +
  ggsci::scale_color_igv(name = 'Phylogenetic order') +
  geom_tippoint(aes(shape = with.rnaseq), color = 'red', size = 5) +
  scale_shape(na.translate = F, name = 'RNA-seq with ≥5 libraries') +
  ggtitle('Chen genomic tree')

################################################################################

p.otu

ggsave('5-rnaseq-overview.jpeg', width = 8, height = 7)


