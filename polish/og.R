
# Script contains adapted plots from
# - B_orthology-groups.R
# - B_seqids.R

library(tidyverse)

###############################################################################
path.tax <- '~/remote-server/projects/rth/co2capture/subprojects/OCyRS/OCyRS-pipeline/data/A_representatives/taxonomy.tsv'
path.genes <- '~/remote-server/projects/rth/co2capture/subprojects/OCyRS/OCyRS-pipeline/data/A_representatives/genes.tsv.gz'

path.og <- '~/remote-server/projects/rth/co2capture/subprojects/OCyRS/OCyRS-pipeline/data/B_OGs.tsv'
path.kegg <- '~/remote-server/projects/rth/co2capture/subprojects/OCyRS/OCyRS-pipeline/data/A_representatives/kegg.tsv.gz'

kegg.url <- 'https://rest.kegg.jp/link/pathway/ko'


################################################################################

tax <- read_tsv(path.tax)
genes <- read_tsv(path.genes)

og <- read_tsv(path.og)
kegg <- read_tsv(path.kegg)

path.ko <- read_tsv(kegg.url, col_names = c('term', 'path')) %>%
  mutate_all(str_remove, '^.*:')


################################################################################
################################################################################
# Inspect phylogeny distribution of order/family/genus/species
# Q: What is the largest sub-clade in which all genomes have at least one
#    copy of the orthologoue

# The genomes per clade
tax %>%
  select(tax.bio, order, family, genus, species) %>%
  gather('sub', 'txid', - tax.bio) %>%
  left_join(
    count(., sub, txid),
    c('sub', 'txid')
  ) %>%
  drop_na -> tax.sub


# check for overlaps
tax.sub %>%
  # filter(n > 3) %>%
  rename(genomes = n) %>%
  arrange(sub, txid) %>%
  left_join(
    kegg %>%
      select(- tax.bio.gene) %>%
      unique,
    'tax.bio'
  ) %>%
  drop_na() %>%
  count(db, term, title, sub, txid, genomes) %>%
  rename(term.genomes = n) -> tax.kegg.lap


# Only keep clades in which all genomes match
tax.kegg.lap %>%
  filter(genomes == term.genomes) %>%
  # largest clade size per term
  group_by(db, term, title) %>%
  summarize(max.clade = max(genomes)) %>%
  ungroup -> kegg.max.clade


################################################################################
# Determine unambigous terms form weighted relative gene count

kegg %>%
  filter(db == 'ko') -> ko

# weight genes by number of different terms they are part of
ko %>%
  count(tax.bio.gene) %>%
  rename(terms = n) %>%
  mutate(w = 1 / terms) %>%
  # sum up weights per term
  left_join(ko, 'tax.bio.gene') %>%
  group_by(term) %>%
  summarize(ws = sum(w)) %>%
  # compare to overall number of genes/geomes of term
  left_join(
    ko %>%
      group_by(term) %>%
      summarize(
        uniq.genes = length(tax.bio.gene),
        uniq.genomes = tax.bio %>%
          unique %>%
          length
      ),
    'term'
  ) %>%
  # add info on what the largest recalled clade was
  left_join(kegg.max.clade, 'term') %>%
  mutate_at('max.clade', replace_na, 0) -> dat.ko

dat.ko %>%
  ggplot(aes(uniq.genomes, uniq.genes)) +
  geom_point() +
  ylab('No. genes annotated in term') +
  xlab('No. of different genomes') +
  theme_bw(18) +
  ggtitle('Sizes of KEGG orthology terms') -> p1

dat.ko %>%
  ggplot(aes(max.clade)) +
  geom_histogram(bins = 60) +
  scale_x_continuous(breaks = seq(0, max(dat.ko$max.clade), 10)) +
  geom_vline(xintercept = c(5, 10, 30), color = 'red') +
  xlab('Largest sub-clade covered by term') +
  theme_bw(18) +
  ggtitle('Overlap of terms with taxonomic sub-clades',
          '(order, family, genus, species)') -> p2

dat.ko %>%
  mutate(
    # ignore cases where NCBI identifies single species and
    # proGenomes has multiple species cluster
    max.clade = ifelse(max.clade == 0, 1, max.clade),
    max.c = cut(max.clade, c(0, 1, 4, 10, 30, Inf),
                include.lowest = TRUE) %>%
      fct_recode(
        '1' = "[0,1]",
        '2..4' = "(1,4]",
        '5..10' = "(4,10]",
        '11.30' = "(10,30]",
        '31+' = "(30,Inf]"
      )
  ) -> foo
foo %>%
  ggplot(aes(uniq.genomes, ws / uniq.genes,
             col = max.c)) +
  scale_color_manual(values = c(
    '#000000', #black
    "#D55E00", #organge
    "#56B4E9", # light blue
    # "#0072B2", # dark blue
    "#009E73", # green
    "#CC79A7"  # pink
  ),
  name = 'Largest sub-clade covered by term'
  ) +
  geom_point(size = 2, alpha = 0.7) +
  xlab('No. of different genomes') +
  ylab('Weighted relative size') +
  theme_bw(18) +
  theme(legend.position = 'bottom') +
  guides(color = guide_legend(nrow = 2)) +
  ggtitle('KEGG orthology term uniqueness') -> p3

foo %>%
  ggplot(aes(max.c, ws / uniq.genes)) +
  geom_violin() +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  geom_hline(yintercept = 0.85, color = 'blue') +
  xlab('Largest sub-clade covered by term') +
  ylab('Weighted relative size') +
  theme_bw(18) +
  ggtitle('Unambiguous orthologies') -> p4


cowplot::plot_grid(
  p1, p2, p3, p4,
  labels = 'AUTO',
  label_size = 18
)
ggsave('supplement-fig_OGs.jpeg', width = 16, height = 10)


################################################################################
# Check for overlaps between pathways and OGs

# 1. the pathway info
kegg %>%
  filter(db == 'pathway') %>%
  select(path = term, pathway = title, tax.bio.gene) -> path
path %>%
  # 2. overlapping OGs
  right_join(og, 'tax.bio.gene') %>%
  left_join(count(og, term), 'term') %>%
  rename(term.size = n) %>%
  # 3. the propotion of OG contained in the pathway
  count(path, pathway, term, title, term.size) %>%
  rename(overlap = n) %>%
  mutate(prop = overlap / term.size) -> path.og

path.og %>%
  drop_na() %>%
  left_join(
    mutate(path.ko, note = 'yes'),
    c('path', 'term')
  ) %>%
  mutate_at('note', replace_na, 'no') %>%
  ggplot(aes(note, prop)) +
  geom_violin(aes(fill = note)) +
  geom_boxplot(width = 0.1) +
  scale_fill_manual(values = c(
    '#56B4E9',
    '#E69F00'
  )) +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  xlab('Curated Pathway-KO association') +
  ylab('Proportion of OG overlapping pathway') +
  theme_bw(18) +
  theme(legend.position = 'hide') -> p1

path.og %>%
  drop_na() %>%
  filter(term.size == overlap) %>%
  select(- overlap, - prop) -> path.og2

path.og2  %>%
  group_by(path, pathway) %>%
  summarize(
    terms = n(),
    avg.ts = mean(term.size)
  ) %>%
  ungroup() %>%
  mutate(
    # highlight pathways of either top5 avg size or number for terms
    l = ifelse(
      (rank(- avg.ts) <= 5) | (rank(- terms) <= 5),
      pathway, NA_character_)
  ) %>%
  ggplot(aes(terms, avg.ts, label = l, color = l)) +
  geom_point(size = 3) +
  ggrepel::geom_label_repel(size = 5) +
  scale_x_log10() +
  xlab('No. OGs in pathway') +
  ylab('Avg. no. genes per OG') +
  theme_bw(18) +
  theme(legend.position = 'hide') -> p2

ggsave(
  plot = cowplot::plot_grid(p1, p2, labels = 'AUTO', label_size = 18),
  filename = 'supplement-fig_path-og.jpeg',
  width = 14, height = 7
)

################################################################################

