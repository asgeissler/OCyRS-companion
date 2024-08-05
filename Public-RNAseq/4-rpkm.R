# Check with the zFPKM method, if any motifs are expressed

library(tidyverse)

library(ggpubr)

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
  write_tsv('4-expression-counts.tsv.gz')

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
  res2 <-
    res |>
    select(- Geneid) |>
    as.data.frame() |>
    magrittr::set_rownames(res$Geneid)
  
  res2
}

expression.rpkm <-
  expression.dat |>
  map(rpkm)



saveRDS(expression.rpkm, '4-rpkm.rds')
