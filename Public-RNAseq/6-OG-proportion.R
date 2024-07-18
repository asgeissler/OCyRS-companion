# Clarify what proporiton of motifs are expressed in comparison to the 
# expression of the ortholog gene anchor

library(tidyverse)

dat <-
  '4-expression-ratios.tsv' |>
  read_tsv() |>
  filter(method == 'RPKM') |>
  select(- method)

kegg <-
  '/home/projects/rth/co2capture/subprojects/OCyRS/OCyRS-pipeline/data/A_representatives/kegg.tsv.gz' |>
  read_tsv()

################################################################################
# match motifs to the gene anchors in the respective genome

motif.gene.pairs <-
  dat |>
  # get motifs
  filter(str_detect(gene, 'fna.motif')) |>
  rename(motif = gene, motif.libs = expressed.libs) |>
  mutate(crs = str_remove(motif, ';.*')) |>
  select(- c(libs.total, ratio)) |>
  # identify tax.bio entry
  separate(genome, c('species', 'txid', 'bio', 'order'), sep = '_') |>
  mutate(tax.bio = paste(
    str_remove(txid, '^txid'),
    bio,
    sep = '.'
  )) |>
  select(- c(txid, bio, order)) |>
  # find the potential ortholog gene anchorS
  # (plural due to paralogs, but a first estimation)
  mutate(term = str_remove(motif, '_.*')) |>
  left_join(
    kegg |>
      filter(db == 'ko') |>
      select(- db),
    c('term', 'tax.bio'),
    relationship = "many-to-many"
  ) |>
  # add expression info of gene
  left_join(
    dat |> select(tax.bio.gene = gene, gene.libs = expressed.libs),
    'tax.bio.gene'
  )

################################################################################
# general stats
  
motif.gene.pairs |>
  select(crs) |>
  unique() |>
  nrow()
# 327 potential CRSs with TA data

motif.gene.pairs |>
  filter(motif.libs >= 3) |>
  select(crs) |>
  unique() |>
  nrow()
# 251 TA in at least one species

################################################################################
# motif in relationship to genes

freqs <-
  motif.gene.pairs |>
  group_by(tax.bio, species, crs) |>
  summarize(
    motif.expressed = max(motif.libs) >= 3,
    gene.expressed = max(gene.libs) >= 3,
  ) |>
  group_by(tax.bio, species) |>
  count(motif.expressed, gene.expressed) |>
  mutate(mode = case_when(
     motif.expressed &  gene.expressed  ~ 'both.expressed',
    !motif.expressed &  gene.expressed  ~ 'gene.expressed',
     motif.expressed & !gene.expressed  ~ 'motif.expressed',
    !motif.expressed & !gene.expressed  ~ 'neither.expressed'
  )) |>
  select(species, mode, n) |>
  spread(mode, n, fill = 0)

# the crs is never expressed if the gene is not
assertthat::assert_that(! 'motif.expressed' %in% colnames(freqs))

# tax.bio             species                                    both.expressed gene.expressed neither.expressed
# 103690.SAMD00061094 Nostoc.sp.PCC.7120.FACHB.418                           20             13                 2
# 1140.SAMN02598254   Synechococcus.elongatus.PCC.7942.FACHB.805             39              2                 0
# 1147.SAMN02471775   Synechocystis.sp.PCC.6714                               4              7                 2
# 1148.SAMD00061113   Synechocystis.sp.PCC.6803                              11              1                 0
# 197221.SAMD00061106 Thermosynechococcus.vestitus.BP.1                       6              3                 0
# 203124.SAMN02598485 Trichodesmium.erythraeum.IMS101                        24              9                 0
# 329726.SAMN02604308 Acaryochloris.marina.MBIC11017                         14             13                 0
# 59920.SAMN00623057  Prochlorococcus.marinus.str.NATL2A                    157              0                 0
# 74546.SAMN02598321  Prochlorococcus.marinus.str.MIT.9312                  136            113                 0

freqs |>
  mutate(
    total = both.expressed + gene.expressed + neither.expressed,
    both.pct = both.expressed / total * 100,
    not.both.pct = (gene.expressed + neither.expressed) / total * 100
  ) |>
# tax.bio             species            both.expressed gene.expressed neither.expressed total both.pct not.both.pct
# 103690.SAMD00061094 Nostoc.sp.PCC.712…             20             13                 2    35     57.1        42.9 
# 1140.SAMN02598254   Synechococcus.elo…             39              2                 0    41     95.1         4.88
# 1147.SAMN02471775   Synechocystis.sp.…              4              7                 2    13     30.8        69.2 
# 1148.SAMD00061113   Synechocystis.sp.…             11              1                 0    12     91.7         8.33
# 197221.SAMD00061106 Thermosynechococc…              6              3                 0     9     66.7        33.3 
# 203124.SAMN02598485 Trichodesmium.ery…             24              9                 0    33     72.7        27.3 
# 329726.SAMN02604308 Acaryochloris.mar…             14             13                 0    27     51.9        48.1 
# 59920.SAMN00623057  Prochlorococcus.m…            157              0                 0   157    100           0   
# 74546.SAMN02598321  Prochlorococcus.m…            136            113                 0   249     54.6        45.4 
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


