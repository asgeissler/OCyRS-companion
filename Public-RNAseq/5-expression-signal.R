# Check with the zFPKM method, if any motifs are expressed

library(tidyverse)
library(patchwork)
library(ggpubr)

library(zFPKM)

################################################################################

expression.rpkm <- readRDS('4-rpkm.rds')

################################################################################
################################################################################
# Compare distributions of RPKM values of genes vs CRS/Rfam


meta <-
  expression.rpkm %>%
  map2(names(.), function(xs, ys) {
    xs |>
      rownames() |>
      tibble(gene = _) |>
      mutate(
        genome = ys,
        type = case_when(
          str_starts(gene, 'RF') ~ 'Rfam',
          str_detect(gene, 'fna.motif')  ~ 'CRS',
          TRUE ~ 'gene'
        )
      )
  }) |>
  bind_rows()

write_tsv(meta, '5-meta.tsv.gz')

################################################################################

meta |>
  count(type, genome) |>
  spread(type, n) |>
  rename(
    CRS.seq = CRS,
    Rfam.hits = Rfam
  ) |>
  left_join(
    meta |>
      filter(type == 'CRS') |>
      mutate(motif = str_remove(gene, ';pos.*')) |>
      select(genome, motif) |>
      unique() |>
      count(genome, name = 'CRS.motifs'),
    'genome'
  ) |> 
  left_join(
    meta |>
      filter(type == 'Rfam') |>
      mutate(fam = str_remove(gene, ';.*')) |>
      select(genome, fam) |>
      unique() |>
      count(genome, name = 'Rfam.fam'),
    'genome'
  ) |> 
  arrange(desc(gene)) |>
  select(genome, gene, Rfam.hits, Rfam.fam, everything()) -> annot.stat
  
annot.stat
# genome                                                                            gene Rfam.hits Rfam.fam CRS.seq CRS.motifs
# <chr>                                                                            <int>     <int>    <int>   <int>      <int>
# 1 Acaryochloris.marina.MBIC11017_txid329726_SAMN02604308_unknown                    6333        90       16      29         27
# 2 Nostoc.sp.PCC.7120.FACHB.418_txid103690_SAMD00061094_Nostocales                   5429       200       28      40         35
# 3 Trichodesmium.erythraeum.IMS101_txid203124_SAMN02598485_Oscillatoriales           4498        58       17      37         33
# 4 Synechocystis.sp.PCC.6714_txid1147_SAMN02471775_Synechococcales                   3565        64       18      13         13
# 5 Synechocystis.sp.PCC.6803_txid1148_SAMD00061113_Synechococcales                   3216        65       20      12         12
# 6 Synechococcus.elongatus.PCC.7942.FACHB.805_txid1140_SAMN02598254_Synechococcales  2663        62       14      41         41
# 7 Thermosynechococcus.vestitus.BP.1_txid197221_SAMD00061106_unknown                 2520        57       16      10          9
# 8 Prochlorococcus.marinus.str.NATL2A_txid59920_SAMN00623057_Synechococcales         2206        55       15     157        157
# 9 Prochlorococcus.marinus.str.MIT.9312_txid74546_SAMN02598321_NA                    2006        58       20     249        249


################################################################################

rdat <-
  expression.rpkm |>
  map(as_tibble, rownames = 'gene') %>%
  map2(names(.), ~ mutate(.x, genome = .y)) |>
  map(pivot_longer, - c(gene, genome)) |>
  bind_rows()

################################################################################
################################################################################

rdat |>
  left_join(meta, c('gene', 'genome')) |>
  # filter(genome == 'Synechocystis.sp.PCC.6714_txid1147_SAMN02471775_Synechococcales') |>
  group_by(genome, name, type) |>
  summarize(
    avg = mean(log10(value + 1)),
    se = sd(log10(value + 1))
  ) |>
  ungroup()  |>
  arrange(genome, type, avg) |>
  mutate_at('name', fct_inorder) |>
  ggplot(aes(
    name, avg,
    ymin = avg - se,
    ymax = avg + se,
    group = type,
    fill = type,
    color = type)) +
  geom_ribbon(alpha = .5) +
  geom_line(size = 1.5) +
  ggsci::scale_color_jama() +
  ggsci::scale_fill_jama() +
  xlab('RNA-seq libraries') +
  ylab('log10(RPKM + 1)') +
  coord_flip() +
  facet_grid(row = 'genome', scales = 'free_y',  space = "free") +
  theme_pubr(18) +
  theme(
    strip.text.y = element_text(angle = 0),
    axis.text.y = element_blank()
  )

ggsave('5-rpkm-avg-se.jpg', width = 20, height = 15)


rdat |>
  left_join(meta, c('gene', 'genome')) |>
  filter(genome == 'Synechocystis.sp.PCC.6714_txid1147_SAMN02471775_Synechococcales') |>
  ggplot(aes(log10(value + 1), color = type)) +
  stat_ecdf() +
  facet_wrap(~ name) +
  theme_pubr(18) +
  ggtitle('Synechocystis.sp.PCC.6714_txid1147_SAMN02471775_Synechococcales')
    
################################################################################
################################################################################
# Check out the zFPKM versions

z.dat <-
  expression.rpkm |>
  map(function(x) {
    x <- as.data.frame(x + 1)
    zFPKM(x)
  })


x <- 'Synechocystis.sp.PCC.6714_txid1147_SAMN02471775_Synechococcales'
z.dat |>
  _[[x]] |>
  as_tibble(rownames = 'gene') %>%
  pivot_longer(- gene) |>
  ggecdf('value', color = 'name', scales = 'free_x') +
  ggsci::scale_color_jco(name = NULL) +
  theme_bw(16) +
  ggtitle(x |> str_remove('_txid.*')) +
  theme(legend.position = 'top') +
  xlab('zFPKM expression') +
  geom_vline(xintercept = -3, color = 'blue') +
  geom_vline(xintercept = -2, color = 'orange') +
  geom_vline(xintercept = -1, color = 'red') +
  scale_y_continuous(breaks = seq(0, 1, .1))

ggsave('5-cutoff.jpeg', width = 18, height = 8)

################################################################################
# save for export to paper

z.dat |>
  map(as_tibble, rownames = 'Geneid') %>%
  map2(names(.), ~ mutate(.x, genome = .y)) |>
  map(pivot_longer, - c(Geneid, genome),
      names_to = 'library', values_to = 'zFPKM') |>
  bind_rows() |>
  dplyr::select(genome, everything()) |>
  write_tsv('5-zFPKM-data.tsv.gz')


################################################################################
# justify cutoff by proportion expressed per library

seq(-3, -0, .25) |>
  map(function(cutoff) {
    z.dat %>%
      map2(names(.), function(x, genome) {
        apply(x > cutoff, 2, sum) %>%
          tibble(
            lib = names(.),
            cutoff = cutoff,
            genome = genome,
            expressed.genes = ., n.genes = dim(x)[1]
          )
      }) |>
          bind_rows()
  }) |>
  bind_rows() -> foo

foo |>
  mutate(prop = expressed.genes / n.genes * 100) |>
  mutate_at('genome', str_remove, '_txid.*') |>
  mutate_at('genome', str_replace_all, '\\.', ' ') |>
  ggplot(aes(cutoff, prop, group = lib, color = genome)) +
  geom_line(alpha = .5) +
  geom_point() +
  scale_color_brewer(palette = 'Paired', type = 'seq') +
  ylab('% Genes detected as expressed\nfor each RNA-seq library') +
  xlab('expression standard deviation below adjusted average (zFPKM threshold)') +
  geom_vline(xintercept = -1, color = 'red') +
  theme_pubr(18) +
  theme(legend.position = 'right')

ggsave('5-per-lib-detected.jpeg', width = 14, height = 8)

################################################################################
# compare choice of normalization on expression detection

z.dat %>%
  map(function(x) {
    apply(x > -1, 1, sum) %>%
      tibble(gene = names(.), expressed.libs = ., libs.total = dim(x)[2]) |>
      mutate(ratio = expressed.libs / libs.total)
  }) %>%
  map2(names(.), ~ mutate(.x, genome = .y)) |>
  bind_rows() -> overall

write_tsv(overall, '5-expression-ratios.tsv')

# overall <- read_tsv('5-expression-ratios.tsv')

overall |>
  mutate_at('genome', str_remove, '_txid.*') |>
  ggplot(aes(ratio, color = genome)) +
  stat_ecdf() +
  xlab('Ratio libraries detected as expressed') +
  guides(color = guide_legend(ncol = 2)) +
  theme_pubr(16)

################################################################################
################################################################################

dat <-
  overall |>
  left_join(meta) |>
  mutate(gene = str_remove(gene, ';.*')) |>
  group_by(gene, genome, type, libs.total) |>
  summarize(max.expressed.libs = max(expressed.libs)) |>
  ungroup() |>
  filter(max.expressed.libs >= 3) |>
  count(genome, type, libs.total) |>
  left_join(
    annot.stat |>
      select(genome, gene, Rfam = Rfam.fam, CRS = CRS.motifs) |>
      pivot_longer(- genome, names_to = 'type')
  ) |>
  mutate(pct = n / value * 100 )

dat |>
  select(genome, type, libs.total, pct) |>
  spread(type, pct)  |>
  rename('Libraries' = libs.total) |>
  GGally::ggpairs(columns = 2:5) +
  theme_bw(18)

ggsave('5-pairs.jpeg', width = 16, height = 12)

################################################################################

genes.total <-
  overall |>
  left_join(meta) |>
  select(gene, genome) |>
  mutate(name = str_remove(Geneid, ';pos[0-9]+$')) |>
  select(genome, name, is.motif) |>
  unique() |>
  group_by(genome) |>
  summarize(
    total.features = n(),
    total.motifs = sum(is.motif)
  )

dat |>
  separate(genome, c('species', 'txid', 'bio', 'order'), sep = '_') |>
  mutate_at('species', str_replace_all, '\\.', ' ') |>
  transmute(
    species, 
    libs.total,
    type,
    total = prettyNum(value, big.mark = ','),
    expressed = sprintf(
      '%s (%s%%)',
      prettyNum(n, big.mark = ','),
      round(pct, 1)
    )
  ) |>
  pivot_longer(c(total, expressed)) |>
  unite('name', c(type, name), sep = '.') |>
  pivot_wider() |>
  select(
    species, libs.total,
    gene.total, gene.expressed,
    Rfam.total, Rfam.expressed,
    CRS.total, CRS.expressed
  ) |>
# species                                    libs.total gene.total gene.expressed Rfam.total Rfam.expressed CRS.total CRS.expressed
# 1 Acaryochloris marina MBIC11017                      9 6,333      5,518 (87.1%)  16         12 (75%)       27        14 (51.9%)   
# 2 Nostoc sp PCC 7120 FACHB 418                       13 5,429      4,966 (91.5%)  28         19 (67.9%)     35        20 (57.1%)   
# 3 Prochlorococcus marinus str MIT 9312                9 2,006      1,912 (95.3%)  20         16 (80%)       249       137 (55%)    
# 4 Prochlorococcus marinus str NATL2A                 73 2,206      2,206 (100%)   15         15 (100%)      157       157 (100%)   
# 5 Synechococcus elongatus PCC 7942 FACHB 805        228 2,663      2,650 (99.5%)  14         11 (78.6%)     41        39 (95.1%)   
# 6 Synechocystis sp PCC 6714                           6 3,565      3,089 (86.6%)  18         14 (77.8%)     13        4 (30.8%)    
# 7 Synechocystis sp PCC 6803                          79 3,216      3,156 (98.1%)  20         18 (90%)       12        11 (91.7%)   
# 8 Thermosynechococcus vestitus BP 1                  26 2,520      2,398 (95.2%)  16         14 (87.5%)     9         6 (66.7%)    
# 9 Trichodesmium erythraeum IMS101                    15 4,498      3,973 (88.3%)  17         14 (82.4%)     33        25 (75.8%)   
  knitr::kable('latex')

################################################################################
################################################################################

overall |>
  filter(expressed.libs >= 3) |>
  left_join(meta) |>
  mutate(gene = str_remove(gene, ';.*')) |>
  select(gene, type, genome) |>
  unique() |>
  dplyr::count(gene, type) |>
  arrange(desc(n)) |>
  dplyr::rename(expressed.in.genomes = n) -> xs

################################################################################

'~/OCyRS-pipeline/data/K_ko-path.tsv' |>
  read_tsv() |>
  select(term, ortholog, pathway) |>
  mutate_at('pathway', replace_na, '') |>
  unique() |>
  group_by(term, ortholog) |>
  summarize(pathways = str_c(pathway, collapse = '; ')) -> xs.meta

xs |>
  mutate(term = str_remove(gene, '_.*')) |>
  inner_join(xs.meta, 'term') -> xs2

write_tsv(xs2, '5-maybe-interest.tsv')

xs2
