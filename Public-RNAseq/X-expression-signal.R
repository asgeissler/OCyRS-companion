# Check with the zFPKM method, if any motifs are expressed

library(tidyverse)
library(patchwork)
library(ggpubr)

library(SummarizedExperiment)
library(zFPKM)

################################################################################
# Double check stats of mapping and assignment rates for sanity check


load.json <- function(i) {
  dat <- jsonlite::read_json(i)
  
  dat$report_general_stats_data %>%
    map(function(xs) {
      map2(xs, names(xs),
           function(x, y) {
             x$sample <- y
             keep(x, ~ length(.x) == 1) %>%
               as_tibble()
           }) %>%
        bind_rows()
    }) %>%
    bind_rows()
}

map.dat <-
  '3-RNA-Browser/*/multiqc/multiqc_data/multiqc_data.json' |>
  Sys.glob() |>
  map(function(i) {
    # i <- '3-RNA-Browser/Acaryochloris.marina.MBIC11017_txid329726_SAMN02604308_unknown/multiqc/multiqc_data/multiqc_data.json'
    j <-
      i |>
      dirname() |>
      dirname() |>
      dirname() |>
      basename()
    i |>
      load.json() |>
      transmute(
        sample,
        'reads after QC' = total_reads,
        'mapping' = total_reads - unpaired_aligned_none,
        'mapping_rate' = mapping / total_reads * 100,
        'mapping_uniquely' = unpaired_aligned_one,
        'mapping_uniquely_rate' = mapping_uniquely / mapping * 100,
        'feature_assigned' = Assigned,
        'feature_assigned_rate' = percent_assigned
      ) |>
      group_by(sample) |>
      summarize_all(~ .x |> discard(is.na) |> unique()) |>
      ungroup() |>
      mutate(genome = j)
  }) |>
  bind_rows()

################################################################################

map.dat |>
  select(genome, sample, contains('rate')) |>
  pivot_longer(contains('rate')) |>
  mutate_at('genome', str_remove, '_txid.*$') |>
  mutate_at('name', fct_inorder) |>
  ggplot(aes(name, value, color = genome)) +
  geom_boxplot() +
  ylab('% reads compared to preceeding step') +
  xlab(NULL) +
  # ggsci::scale_color_jco(name = NULL) +
  scale_color_viridis_d(name = NULL) +
  guides(color = guide_legend(ncol = 2)) +
  theme_pubr(18)

ggsave('4-mapping-stat.jpeg', width = 12, height = 8)

################################################################################
# Load expression data

expression.dat <-
  '3-RNA-Browser/*/analysis/41_counts.tsv' |>
  Sys.glob() %>%
  set_names(
    . |>
      dirname() |>
      dirname() |>
      basename()
  ) |>
  map(read_tsv)

################################################################################
# save for export to paper

expression.dat %>%
  map2(names(.), ~ mutate(.x, genome = .y)) |>
  map(pivot_longer, -c(genome, Geneid, Chr, Start, End, Strand, Length),
      names_to = 'library', values_to = 'expression_count') |>
  bind_rows() |>
  dplyr::select(genome, everything()) |>
  write_tsv('4-expression-data.tsv.gz')

################################################################################
# prepare FPKM normalized matrices


rpkm <- function(x) {
  # x <- expression.dat$Synechocystis.sp.PCC.6803_txid1148_SAMD00061113_Synechococcales
  x.long <-
    x |>
    select(- c(Chr, Start, End, Strand)) |>
    pivot_longer(- c(Geneid, Length))
  
  x.long |>
    left_join(
      map.dat |>
        select(name = sample, total = mapping_uniquely),
      # x.long |>
      #   group_by(name) |>
      #   summarize(total = sum(value)),
      'name'
    ) |>
    mutate(
      pm = total / 1e6,
      klen = Length / 1e3,
      rpkm = value / pm / klen
    ) |>
    select(Geneid, name, value = rpkm) |>
    pivot_wider() -> res
  res <-
    res |>
    select(- Geneid) |>
    as.data.frame() |>
    magrittr::set_rownames(res$Geneid)
  
  # mask <- rowMax(as.matrix(res)) > 0
  # res[mask, ]
  res
}

expression.rpkm <-
  expression.dat |>
  map(rpkm)

################################################################################
# prepare TPM normalized matrices


tpm <- function(x) {
  # x <- expression.dat$Synechocystis.sp.PCC.6803_txid1148_SAMD00061113_Synechococcales
  x.long <-
    x |>
    select(- c(Chr, Start, End, Strand)) |>
    pivot_longer(- c(Geneid, Length))
  
  step1 <-
    x.long |>
    mutate(
      klen = Length / 1e3,
      v2 = value / klen
    )
  step1 |>
    left_join(
      step1 |>
        group_by(name) |>
        summarise(pm = sum(v2) / 1e6),
      'name'
    ) |>
    transmute(
      Geneid, name,
      value = v2 / pm
    ) |>
    pivot_wider() -> res
  res <-
    res |>
    select(- Geneid) |>
    as.data.frame() |>
    magrittr::set_rownames(res$Geneid)
  
  # mask <- rowMax(as.matrix(res)) > 0
  # res[mask, ]
  res
}

expression.tpm <-
  expression.dat |>
  map(tpm)

################################################################################

my.des <- function(x) {
  # x <- expression.dat$Synechocystis.sp.PCC.6803_txid1148_SAMD00061113_Synechococcales
  x.mat <-
    x |>
    select(- c(Geneid, Chr, Start, End, Strand, Length)) |>
    as.matrix() |>
    magrittr::set_rownames(x$Geneid)
  
  # mask <- rowMax(x.mat) > 0
  
  DESeq2::DESeqDataSetFromMatrix(
    # x.mat[mask, ],
    x.mat,
    tibble(lib = colnames(x.mat)),
    ~ 1
  ) |>
    DESeq2::DESeq()
}

des.list <-
  expression.dat |>
  map(my.des)
  
norm.versions <-
  list(
    RPKM = expression.rpkm,
    TPM = expression.tpm,
    VST = des.list |>
      map(DESeq2::vst) |>
      map(SummarizedExperiment::assay) |>
      map(as.data.frame) |>
      # avoid negative values in vst
      map(~ .x - min(.x))
  )

################################################################################

norm.versions$`normalized counts` <-
  des.list |>
  map(DESeq2::counts, normalized = TRUE)

all.lengths <-
  expression.dat |>
  map(select, Geneid, Length) |>
  bind_rows() |>
  with(set_names(Length, Geneid))
norm.versions$`normalized counts / gene length` <-
  des.list |>
  map(DESeq2::counts, normalized = TRUE) |>
  map(~ .x / all.lengths[rownames(.x)])
# matrix(1:9, nrow = 3) / 1:3

# norm.versions$`VST / gene length` <-
#   norm.versions$VST |>
#   map(~ .x / all.lengths[rownames(.x)])
    
################################################################################
################################################################################
# Justify choice of normalization method

x <- 'Synechocystis.sp.PCC.6714_txid1147_SAMN02471775_Synechococcales'

norm.versions |>
  map(x) |>
  map(as_tibble, rownames = 'gene') %>%
  map2(names(.), ~ mutate(.x, method = .y)) |>
  bind_rows() |>
  pivot_longer(- c(gene, method)) -> foo

foo |>
  mutate_at('value', ~ .x + 1) |>
  ggplot(aes(value, color = name)) +
  geom_density() +
  scale_x_log10() +
  facet_wrap(~ method, scales = 'free') +
  ggsci::scale_color_d3(name = NULL) +
  xlab('Expression value after tranformation') +
  theme_pubr(18) +
  ggtitle(x |> str_remove('_txid.*'))

ggsave('4-norm-choice.jpeg', width = 18, height = 8)

norm.versions |>
  map(x) %>%
  map2(names(.), ~ tibble(
    value = apply(.x, 1, median),
    gene = rownames(.x),
    length = all.lengths[gene],
    method = .y
  )) |>
  bind_rows() -> bar
bar |>
  mutate_at('value', ~ .x + 1) |>
  ggscatter(
    'length', 'value',
    facet.by = 'method',
    scales =  'free',
    cor.coef = TRUE,
    cor.coeff.args = list(color = 'red', size = 5),
    add = 'reg.line',
    add.params = list(color = 'red')
  ) +
  scale_x_log10() +
  scale_y_log10() +
  ggtitle(x |> str_remove('_txid.*')) +
  xlab('Gene length') +
  ylab('Gene expression')

ggsave('4-norm-length-correlation.jpeg', width = 14, height = 8)

################################################################################


################################################################################
# Check out the various zFPKM data versions

z.dat <-
  norm.versions %>%
  map(function(mats) {
    mats |>
      map(function(x) {
        x <- as.data.frame(x + 1)
        zFPKM(x)
      })
  })


x <- 'Synechocystis.sp.PCC.6714_txid1147_SAMN02471775_Synechococcales'
z.dat |>
  map(x) |>
  map(as_tibble, rownames = 'gene') %>%
  map2(names(.), ~ mutate(.x, norm = .y)) |>
  bind_rows() |>
  pivot_longer(- c( gene, norm)) |>
  ggecdf('value', color = 'name', facet.by = 'norm', scales = 'free_x') +
  ggsci::scale_color_jco(name = NULL) +
  theme_bw(16) +
  ggtitle(x |> str_remove('_txid.*')) +
  theme(legend.position = 'top') +
  xlab('zFPKM expression') +
  geom_vline(xintercept = -3, color = 'blue') +
  geom_vline(xintercept = -2, color = 'orange') +
  geom_vline(xintercept = -1, color = 'red') +
  # scale_x_continuous(breaks = -4:4) +
  scale_y_continuous(breaks = seq(0, 1, .1))

ggsave('4-new-cutoff.jpeg', width = 18, height = 8)

################################################################################
# save for export to paper

z.dat$RPKM |>
  map(as_tibble, rownames = 'Geneid') %>%
  map2(names(.), ~ mutate(.x, genome = .y)) |>
  map(pivot_longer, - c(Geneid, genome),
      names_to = 'library', values_to = 'zFPKM') |>
  bind_rows() |>
  dplyr::select(genome, everything()) |>
  write_tsv('4-zFPKM-data.tsv.gz')


################################################################################
# justify cutoff by proportion expressed per library

seq(-3, -0, .25) |>
  map(function(cutoff) {
    z.dat %>%
      map2(names(.), function(mats, method) {
        mats %>%
          map2(names(.), function(x, genome) {
            apply(x > cutoff, 2, sum) %>%
              tibble(
                lib = names(.),
                cutoff = cutoff,
                method = method,
                genome = genome,
                expressed.genes = ., n.genes = dim(x)[1]
              )
          }) |>
          bind_rows()
    }) |>
      bind_rows()
  }) |>
  bind_rows() -> foo

foo |>
  filter(method == 'RPKM') |>
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
  # facet_wrap(~ method)

ggsave('4-per-lib-detected.jpeg', width = 14, height = 8)

################################################################################
# compare choice of normalization on expression detection

z.dat %>%
  map2(names(.), function(mats, method) {
    mats |>
      map(function(x) {
        apply(x > -1, 1, sum) %>%
          tibble(gene = names(.), expressed.libs = ., libs.total = dim(x)[2]) |>
          mutate(ratio = expressed.libs / libs.total)
      }) %>%
      map2(names(.), ~ mutate(.x, genome = .y)) |>
      bind_rows() |>
      mutate(method = method)
  }) |>
  bind_rows() -> overall

write_tsv(overall, '4-expression-ratios.tsv')

# overall <- read_tsv('4-expression-ratios.tsv')

overall |>
  mutate_at('genome', str_remove, '_txid.*') |>
  ggplot(aes(ratio, color = genome)) +
  stat_ecdf() +
  facet_wrap(~ method) +
  xlab('Ratio libraries detected as expressed') +
  guides(color = guide_legend(ncol = 2)) +
  theme_pubr(16)
  # theme(legend.position = 'right')

################################################################################

overall |>
  # filter(ratio >= 0.5) |>
  filter(expressed.libs >= 3) |>
  select(method, genome, gene) %>%
  split(.$genome) |>
  map(~ split(.x$gene, .x$method)) -> foo

overall %>%
  split(.$genome) |>
  map(function(x) { foo[[dplyr::first(x$genome)]]$Genome <<- x$gene})
  

foo %>%
  map2(names(.), function(xs, genome) {
    venn::venn(
      xs[c('Genome', 'RPKM', 'TPM', 'normalized counts')],
      ilabels = 'counts',
      ilcs = 1.3, sncs = 1.3,
      box = FALSE,
      zcolor = 'style',
      ggplot = TRUE
    ) +
      annotate('text', 50, 980,
               label = genome |> str_remove('_txid.*'),
               hjust = 0, size = 5)
  }) |>
  purrr::reduce(.f = `+`)

ggsave('4-methods-venn.jpeg', width = 20, height = 20)


################################################################################

genes.total <-
  expression.dat |>
  map(select, Geneid) %>%
  map2(names(.), ~ mutate(.x, genome = .y)) |>
  bind_rows() |>
  mutate(is.motif = str_detect(Geneid, 'fna.motif')) |>
  mutate(name = str_remove(Geneid, ';pos[0-9]+$')) |>
  select(genome, name, is.motif) |>
  unique() |>
  group_by(genome) |>
  summarize(
    total.features = n(),
    total.motifs = sum(is.motif)
  )

overall |>
  filter(method == 'RPKM') |>
  filter(expressed.libs >= 3) |>
  mutate(is.motif = str_detect(gene, 'fna.motif')) |>
  mutate(name = str_remove(gene, ';pos[0-9]+$')) |>
  select(name, is.motif, libs.total, genome) |>
  unique() |>
  group_by(genome) |>
  summarize(
    expressed.features = n(),
    expressed.motifs = sum(is.motif),
    libs.total = unique(libs.total)
  ) |>
  left_join(genes.total, 'genome') |>
  separate(genome, c('species', 'txid', 'bio', 'order'), sep = '_') |>
  mutate_at('species', str_replace_all, '\\.', ' ') |>
  mutate(
    total.genes = total.features - total.motifs,
    expressed.genes = expressed.features - expressed.motifs
  ) |>
  transmute(
    species, 
    libs.total,
    totalgenes = prettyNum(total.genes, big.mark = ','),
    expressed.genes = sprintf(
      '%s (%s%%)',
      prettyNum(expressed.genes, big.mark = ','),
      round(expressed.genes / total.genes * 100, 1)
    ),
    total.motifs,
    expressed.motifs = sprintf(
      '%s (%s%%)',
      prettyNum(expressed.motifs, big.mark = ','),
      round(expressed.motifs / total.motifs * 100, 1)
    )
  ) |>
# species                                    libs.total totalgenes expressed.genes total.motifs expressed.motifs
# 1 Acaryochloris marina MBIC11017                      9 6,333      5,513 (87.1%)             27 14 (51.9%)      
# 2 Nostoc sp PCC 7120 FACHB 418                       13 5,429      4,960 (91.4%)             35 20 (57.1%)      
# 3 Prochlorococcus marinus str MIT 9312                9 2,006      1,904 (94.9%)            249 136 (54.6%)     
# 4 Prochlorococcus marinus str NATL2A                 73 2,206      2,206 (100%)             157 157 (100%)      
# 5 Synechococcus elongatus PCC 7942 FACHB 805        228 2,663      2,650 (99.5%)             41 39 (95.1%)      
# 6 Synechocystis sp PCC 6714                           6 3,565      3,073 (86.2%)             13 4 (30.8%)       
# 7 Synechocystis sp PCC 6803                          79 3,216      3,156 (98.1%)             12 11 (91.7%)      
# 8 Thermosynechococcus vestitus BP 1                  26 2,520      2,382 (94.5%)              9 6 (66.7%)       
# 9 Trichodesmium erythraeum IMS101                    15 4,498      3,974 (88.4%)             33 24 (72.7%)    
  knitr::kable('latex')

################################################################################
################################################################################

overall |>
  # filter(method == 'normalized counts / gene length') |>
  filter(method == 'RPKM') |>
  # filter(ratio >= .5) |>
  filter(expressed.libs >= 3) |>
  mutate(is.motif = str_detect(gene, 'fna.motif')) |>
  group_by(genome) |>
  summarize(
    expressed.genes = n(),
    expressed.motifs = sum(is.motif)
  ) |>
  left_join(genes.total, 'genome') |>
  mutate(
    total.genes = total.features - total.motifs,
    frac.genes = expressed.genes / total.genes,
    frac.motifs = expressed.motifs / total.motifs,
    lab = sprintf(
      '%s\n%s of %s motif positions',
      str_remove(genome, '_txid.*'),
      expressed.motifs, total.motifs
    )
  ) |>
  ggplot(aes(frac.genes, frac.motifs, color = genome, label = lab)) +
  geom_point() +
  ggrepel::geom_text_repel(size = 5) +
  xlab('Fraction of genes detected as expressed') +
  ylab('Fraction of motifs detected as expressed') +
  theme_pubr(18) +
  theme(legend.position = 'hide') -> p1

overall |>
  # filter(method == 'normalized counts / gene length') |>
  filter(method == 'RPKM') |>
  # filter(ratio >= .5) |>
  filter(expressed.libs >= 3) |>
  mutate(is.motif = str_detect(gene, 'fna.motif')) |>
  group_by(genome, libs.total) |>
  summarize(
    expressed.genes = n(),
    expressed.motifs = sum(is.motif)
  ) |>
  left_join(genes.total, 'genome') |>
  mutate(
    total.genes = total.features - total.motifs,
    frac.genes = expressed.genes / total.genes,
    frac.motifs = expressed.motifs / total.motifs,
    lab = sprintf(
      '%s\n%s of %s motif positions',
      str_remove(genome, '_txid.*'),
      expressed.motifs, total.motifs
    )
  ) |>
  ggplot(aes(frac.genes, libs.total, color = genome, label = lab)) +
  geom_point() +
  ggrepel::geom_text_repel(size = 5) +
  scale_y_log10() +
  annotation_logticks(sides = 'l') +
  xlab('Fraction of genes detected as expressed') +
  ylab('# RNA-seq libraries') +
  theme_pubr(18) +
  theme(legend.position = 'hide') -> p2


p1 + p2
    
ggsave('4-expressed.jpeg', width = 20, height = 10)
           
    
################################################################################

overall |>
  # filter(method == 'normalized counts / gene length') |>
  filter(method == 'RPKM') |>
  # filter(ratio >= .5) |>
  filter(expressed.libs >= 3) |>
  mutate(is.motif = str_detect(gene, 'fna.motif')) |>
  filter(is.motif) |>
  mutate(motif = str_remove(gene, ';pos.*')) |>
  select(motif, genome) |>
  unique() |>
  dplyr::count(motif) |>
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
  mutate(term = str_remove(motif, '_.*')) |>
  left_join(xs.meta, 'term') -> xs2
View(xs2)

write_tsv(xs2, '4-maybe-interest.tsv')

xs2
