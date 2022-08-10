# Query entrez for public RNA-seq data

library(tidyverse)
library(rentrez)
library(xml2)

library(ape)
library(ggtree)

library(patchwork)
library(ggpubr)

motif.tax <-
  '~/OCyRS-pipeline/data/K_motif-tax.tsv' |>
  read_tsv()

tax.n <-
  motif.tax |>
  select(motif, tax.bio) |>
  separate(tax.bio, c('tax', 'bio'), sep = '\\.') |>
  select(motif, tax) |>
  unique() |>
  count(tax) |>
  arrange(desc(n))

################################################################################

# Search for single-end Illumina RNA-seq for all taxonomies that contain
# predicted RNA structure motifs
query <- 'txid%s[Organism] AND ("biomol rna"[Properties] AND "platform illumina"[Properties]) AND "library layout single"[Properties]'
bio.ids <-
  tax.n |>
  pull(tax) |>
  unique() %>%
  set_names(.) |>
  map(sprintf, fmt = query) |>
  map(entrez_search, db = 'sra', retmax = 900)

# retrieve summaries for all datasets
# with extra info, eg no. libraries
bio.summaries <-
  bio.ids |>
  # exclude biosamples without ids
  keep(~ length(.x$ids) > 0) |>
  # query
  map(~ entrez_summary(.x$ids, db = 'sra', always_return_list = TRUE))

# parse entrez summary as tibble
helper.xml <- function(i) {
  # i <- bio.summaries$`1219`
  xml_dat <- i$expxml |>
    sprintf(fmt = '<root>%s</root>') |>
    read_xml()
  
  # Custom parse entrez
  res <-
    tribble(
      ~ name, ~ path, ~ attr,
      'title', '/root/Summary/Title', 'RAW_TEXT',
      'instrument', '/root/Summary/Platform', 'instrument_model',
      'total.spots', '/root/Summary/Statistics', 'total_spots',
      'total.bases', '/root/Summary/Statistics', 'total_bases',
      'submitter', '/root/Submitter', 'center_name',
      'SRA', '/root/Submitter', 'acc',
      'SRP', '/root/Study', 'acc',
      'tax', '/root/Organism', 'taxid',
      'tax.name', '/root/Organism', 'ScientificName',
      'strategy', '/root/Library_descriptor/LIBRARY_STRATEGY', 'RAW_TEXT',
      'is.single', '/root/Library_descriptor/LIBRARY_LAYOUT/SINGLE', 'RAW_EXISTS',
      'is.paired', '/root/Library_descriptor/LIBRARY_LAYOUT/PAIRED', 'RAW_EXISTS'
    ) |>
    pmap(function(name, path, attr) {
      # path <- '/root/Summary/Title'
      node <- xml_find_first(xml_dat, path)
      case_when(
        attr == 'RAW_TEXT' ~ xml_text(node),
        attr == 'RAW_EXISTS' ~ as.character(length(node) > 0),
        TRUE ~ xml_attr(node, attr)
      ) |>
        set_names(name)
    }) |>
    unlist()
  # add run identifies
  res$SRR <-
    i$runs |>
    str_extract_all('(?<=Run acc=")[A-Z]+[0-9]+') |>
    unlist() |>
    str_c(collapse = ';')
  return(res)
}

length(bio.summaries)

rnaseq <-
  bio.summaries |>
  map(~ map(.x, helper.xml)) |>
  map(bind_rows) %>%
  map2(names(.), ~ mutate(.x, tax.search = .y)) |>
  bind_rows()

  
# consistent with query all are single end
assertthat::assert_that(all( rnaseq |> pull(is.single) ))

write_tsv(rnaseq, '0-public-rnaseq-projects.tsv')
# rnaseq <- read_tsv('0-public-rnaseq-projects.tsv') |> mutate_at('tax', as.character)

################################################################################

rnaseq |>
  count(tax == tax.search)
 # FALSE                 101
 # TRUE                  927
# cases other strain than what was searched for -> ignore
rna2 <- 
  rnaseq |>
  filter(tax == tax.search)


# no. libraries etc
stat.overview <-
  full_join(
    rna2 |>
      count(tax, name = 'no.libs'),
    rna2 |>
      select(tax, SRP) |>
      unique() |>
      count(tax, name = 'no.datasets'),
    'tax'
  ) |>
  left_join(tax.n, 'tax') |>
  # add name
  left_join(
    rna2 |>
      select(tax, tax.name) |>
      drop_na() |>
      unique() |>
      group_by(tax) |>
      summarize(species = str_c(tax.name, collapse = '\n')),
    'tax'
  )

p.rnaseq <-
  stat.overview |>
  mutate(
    avg = no.libs / no.datasets,
    nice = ifelse((no.libs >= 0) & (no.datasets >= 2), species, NA_character_)
  ) |>
  ggplot(aes(n, no.libs, color = as.factor(no.datasets), label = nice)) +
  geom_point(size = 5) +
  scale_color_viridis_d(name = 'No. SRA projects (datasets)') +
  scale_y_log10() +
  annotation_logticks(sides = 'l') +
  geom_hline(yintercept = 5, color = 'red') +
  # geom_vline(xintercept = 100, color = 'red') +
  ggrepel::geom_label_repel(color = 'black', nudge_y = .5, force = 5, show.legend = FALSE) +
  xlab('No. candidate RNA structure motifs in genome') +
  ylab('No. RNA-seq libraries, sum across all datasets') +
  ggtitle('Public single-end Illumina RNA-seq data', 'entrez queries on txid') +
  theme_pubr(18)



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
  

x.tax$`RNAseq with ≥5 libs` <-
  stat.overview |>
  filter(no.libs >= 5) |>
  pull(tax)


p.venn <- venn::venn(x.tax, zcolor = 'style', ilcs = 1, sncs = 1.2, ggplot =  TRUE)

################################################################################
# location of species with rnaseq for selection in species tree

ref <- ref.trees$Chen_genomic
ref$tip.label <- str_remove(ref$tip.label, '\\..*$')

# ncbi order for coloration to accept 
rec.orders <-
  '~/OCyRS-pipeline/data/A_representatives/taxonomy.tsv' |>
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
  rename(order = ranks) |>
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
  stat.overview |>
  filter(no.libs >= 5) |>
  mutate(with.rnaseq = 'With RNA-seq ≥5 libs')

p.tree <- ggtree(ref, layout = 'circular', branch.length = 'none') %<+% rna.extra


p.otu <-
  groupOTU(p.tree, ncbi.orders, 'Order') +
  aes(color = Order) +
  ggsci::scale_color_igv(name = 'Phylogenetic order') +
  geom_tippoint(aes(shape = with.rnaseq), color = 'red', size = 5) +
  scale_shape(na.translate = F, name = NULL) +
  ggtitle('Chen genomic tree')

################################################################################


p.full <- (p.rnaseq | (wrap_elements(full = p.venn) / wrap_elements(full = p.otu)))

ggsave('0-overview.jpeg', p.full, width = 16, height = 14)


################################################################################
# save stats and orders for subsequent processing

tbl.orders <-
  ncbi.orders |>
  map(~ tibble(tax = .x)) %>%
  map2(names(.), ~ mutate(.x, order = .y)) |>
  bind_rows()

stat.overview |>
  left_join(tbl.orders, 'tax') |>
  write_tsv('0-public-stat.tsv')

