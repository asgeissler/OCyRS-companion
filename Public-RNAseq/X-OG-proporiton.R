# Clarify what proporiton of motifs are expressed in comparison to the 
# expression of the ortholog gene anchor

library(tidyverse)

library(plyranges)

################################################################################

dat <-
  '4-expression-ratios.tsv' |>
  read_tsv() |>
  filter(method == 'RPKM') |>
  select(- method)

motifpos <-
  '/home/projects/rth/co2capture/subprojects/OCyRS/OCyRS-pipeline/data/J_motif-aln-seq-pos.tsv' |>
  read_tsv()
  

regions <-
  '/home/projects/rth/co2capture/subprojects/OCyRS/OCyRS-pipeline/data/J_novel/all_intergenic_regions.tsv.gz' |>
  read_tsv()

kegg <-
  '/home/projects/rth/co2capture/subprojects/OCyRS/OCyRS-pipeline/data/A_representatives/kegg.tsv.gz' |>
  read_tsv()

novel.crs <-
  '/home/projects/rth/co2capture/subprojects/OCyRS/OCyRS-pipeline/data/J_novel/potentially-novel-motifs.tsv' |>
  read_tsv()

novel.red <-
  '/home/projects/rth/co2capture/subprojects/OCyRS/OCyRS-pipeline/data/L_redundant.tsv' |>
  read_tsv()

################################################################################
# match motifs to the gene anchors in the respective genome


mposrange <-
  motifpos |>
  transmute(
    seqnames = chr,
    start, end, strand,
    motif,
    mid = sprintf('%s;pos%s', motif, start)
  ) |>
  as_granges()

regionsrange <-
  regions |>
  as_granges()


mr.match <-
  join_overlap_intersect(
    mposrange, regionsrange
  ) |>
  mutate(mregion = str_remove(mid, '.fna.*$')) |>
  filter(mregion == region) |>
  as_tibble() |>
  select(motif, mid, origin.gene) |>
  unique()

################################################################################
# create overview stat table with on/off expressed libs

motif.gene.pairs <-
  mr.match |>
  inner_join(dat, c('mid' = 'gene')) |>
  left_join(dat, c('origin.gene' = 'gene')) |>
  select(
    motif, mid, motif.libs = expressed.libs.x,
    gene = origin.gene, gene.libs = expressed.libs.y
  ) |>
  mutate_at('gene.libs', replace_na, 0)


################################################################################
# general stats
  
motif.gene.pairs |>
  select(motif) |>
  unique() |>
  nrow()
# 327 potential CRSs with TA data

motif.gene.pairs |>
  filter(motif.libs >= 3) |>
  select(motif) |>
  unique() |>
  nrow()
# 251 TA in at least one species

################################################################################

motif.gene.pairs |>
  mutate(
    motif.expressed = motif.libs >= 3,
    gene.expressed  = gene.libs >= 3,
    mode = case_when(
       motif.expressed &  gene.expressed  ~ 'both.expressed',
      !motif.expressed &  gene.expressed  ~ 'gene.expressed',
       motif.expressed & !gene.expressed  ~ 'motif.expressed',
      !motif.expressed & !gene.expressed  ~ 'neither.expressed'
    )
  ) -> foo

foo |>
  count(motif, name = 'pos') |>
  left_join(
    foo |>
      filter(mode == 'both.expressed') |>
      count(motif),
    'motif'
  ) |>
  mutate(
    n = replace_na(n, 0),
    ratio = n / pos
  ) |>
  filter(ratio == 1) |>
  nrow()

################################################################################
# motif in relationship to genes

motif.gene.pairs |>
  arrange(mid) |>
  View()

freqs <-
  motif.gene.pairs |>
  mutate(txid = str_remove(gene, '\\..*$')) |>
  group_by(motif, txid) |>
  mutate(
    motif.expressed = max(motif.libs) >= 3,
    gene.expressed  = max(gene.libs) >= 3,
  ) |>
  mutate(mode = case_when(
     motif.expressed &  gene.expressed  ~ 'both.expressed',
    !motif.expressed &  gene.expressed  ~ 'gene.expressed',
     motif.expressed & !gene.expressed  ~ 'motif.expressed',
    !motif.expressed & !gene.expressed  ~ 'neither.expressed'
  )) |>
  ungroup()


freqs |>
  anti_join(freqs |> filter(mode != 'both.expressed'), 'motif') |>
  select(motif) |> 
  unique() |>
  nrow()




# the crs is never expressed if the gene is not
# assertthat::assert_that(! 'motif.expressed' %in% colnames(freqs))
# NOT given


freqs |>
  mutate(
    total = both.expressed + gene.expressed + motif.expressed  + neither.expressed,
    both.pct = both.expressed / total * 100,
    not.both.pct = (gene.expressed + neither.expressed) / total * 100
  ) |>
  View()
  pull(not.both.pct) |>
  # summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   8.333  33.333  31.048  45.382  69.231 
  sd()
# 23.12877

################################################################################

dat |>
  count(genome, libs.total, name = 'genes.total')  |>
  left_join(
    dat |>
      filter(expressed.libs >= 3) |>
      count(genome, name = 'genes.expressed'),
    'genome'
  ) |>
  mutate(prop = genes.expressed / genes.total) |>
  ggpubr::ggscatter('prop', 'libs.total', add = 'reg.line', cor.coef = TRUE)



################################################################################
################################################################################
# Ongoing extra requests
  
mr.match |>
  semi_join(novel.crs, 'motif') |>
  anti_join(novel.red, c('motif' = 'redundant')) |>
  inner_join(dat, c('mid' = 'gene')) |>
  mutate(
    reg = str_remove(motif, '.fna.*$'),
    txid = str_remove(origin.gene, '\\.[^.]*$')
  ) -> foo

foo |>
  filter(expressed.libs >= 3) |>
  # count(motif, txid, reg) |>
  count(motif, txid, origin.gene) |>
  filter(n > 1) |>
  pull(motif) |> unique() |> length()

################################################################################

mr.match |>
  semi_join(novel.crs, 'motif') |>
  anti_join(novel.red, c('motif' = 'redundant')) |>
  inner_join(dat, c('mid' = 'gene')) |>
  mutate(
    reg = str_remove(motif, '.fna.*$'),
    txid = str_remove(origin.gene, '\\.[^.]*$')
  ) -> foo

foo |>
  filter(expressed.libs >= 3) |>
  select(motif, txid, origin.gene) |>
  unique() |>
  count(motif, txid) |>
  filter(n > 1) |>
  pull(motif) |> unique() |> length()


################################################################################

mr.match |>
  semi_join(novel.crs, 'motif') |>
  anti_join(novel.red, c('motif' = 'redundant')) |>
  mutate(
    og = str_remove(motif, '_.*$'),
    genome = str_remove(origin.gene, '\\.[^.]*$')
  ) |>
  select(- mid) |>
  unique() -> foo

foo |>
  count(motif, og, genome) |>
  filter(n > 1) |>
  # select(og) |> unique() |> nrow()
  # 45
  select(motif) |>
  unique() |>
  # nrow()
  # 56
  left_join(foo, 'motif')  |>
  left_join(kegg, c('origin.gene' = 'tax.bio.gene', 'genome' = 'tax.bio'),
            relationship = "many-to-many") -> baz

baz |>
  select(motif, origin.gene, genome) |>
  unique() |>
  count(motif, genome, name = 'og.genes.in.genome') |>
  filter(og.genes.in.genome > 1) -> bar

baz  |>
  filter(og != term) |>
  count(motif, genome, db, term, title, name = 'term.freq') |>
  inner_join(bar, c('motif', 'genome')) |>
  mutate(ratio = term.freq / og.genes.in.genome) |>
  arrange(desc(ratio)) -> fooo

fooo |>
  filter(db == 'module') |>
  filter(term.freq >= 3) |>
  select(motif, genome, db, term, title) |>
  unique() |>
  left_join(baz) |>
  left_join(kegg, c('origin.gene' = 'tax.bio.gene')) |>
  filter(db.y != 'pathway') |>
  View()
