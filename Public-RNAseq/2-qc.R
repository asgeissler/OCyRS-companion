# 1. Process RNA-Schlange results
# 2. QC filter libraries
# 3. Select datasets for subsequent RNA-browser bowtie2 mapping aso

library(tidyverse)
library(ggpubr)

library(patchwork)

################################################################################
# helper to read MultiQC json

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

qc.dat <-
  '1-RNA-Schlange/*/analysis/53_main_multiqc/multiqc_data/multiqc_data.json' |>
  Sys.glob() |>
  map(function(i) {
    # i <- '1-RNA-Schlange/Pleurocapsa.sp.PCC.7319_txid118161_SAMN02261343_Pleurocapsales/analysis/53_main_multiqc/multiqc_data/multiqc_data.json'
    j <-
      i |>
      dirname() |>
      dirname() |>
      dirname() |>
      dirname() |>
      basename()
    i |>
      load.json() |>
      select(
        sample,
        before_filtering_total_reads,
        after_filtering_total_reads,
        non_rRNA,
        num_mapped,
        library_types
      ) |>
      unnest(library_types, keep_empty = TRUE) |>
      group_by(sample) |>
      summarize_all(~ .x |> discard(is.na) |> unique()) |>
      ungroup() |>
      mutate(genome = j)
  }) |>
  bind_rows()

################################################################################
# Salmon predicted library types
qc.dat |>
  count(library_types) |>
  mutate(x = n / sum(n) * 100)

write_tsv(qc.dat, '2-RNA-Schlange.tsv')

################################################################################

p1 <-
  qc.dat |>
  mutate_at('genome', str_replace_all, '_', '\n') |>
  mutate_at('genome', str_replace_all, '\nSAM', ' SAM') |>
  mutate_at('genome', ~ paste(.x, 'order')) |>
  ggplot(aes(genome, fill = library_types)) +
  geom_bar(position = 'stack')  +
  ggsci::scale_fill_jama() +
  xlab(NULL) +
  ylab('RNA-seq libraries') +
  coord_flip() +
  theme_pubr(18)

################################################################################

p2 <-
  qc.dat |>
  mutate_at('genome', str_replace_all, '_', '\n') |>
  mutate_at('genome', str_replace_all, '\nSAM', ' SAM') |>
  mutate_at('genome', ~ paste(.x, 'order')) |>
  ggplot(aes(genome, before_filtering_total_reads)) +
  geom_boxplot() +
  scale_y_log10() +
  xlab(NULL) +
  ylab('No. reads') +
  coord_flip() +
  theme_pubr(18)

################################################################################

p3 <-
  qc.dat |>
  transmute(
    sample, genome,
    filtering = after_filtering_total_reads  / before_filtering_total_reads,
    mRNA = non_rRNA / after_filtering_total_reads,
    salmon_mapped = num_mapped / non_rRNA
  ) |>
  mutate_if(is.numeric, ~ .x * 100) |>
  pivot_longer(- c(genome, sample)) |>
  mutate_at('name', fct_relevel, 'filtering', 'mRNA', 'salmon_mapped') |>
  mutate_at('genome', str_replace_all, '_', '\n') |>
  mutate_at('genome', str_replace_all, '\nSAM', ' SAM') |>
  mutate_at('genome', ~ paste(.x, 'order')) |>
  ggplot(aes(name, value, color = genome)) +
  geom_boxplot() +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  ylab('% reads compared to preceeding step') +
  xlab('QC processing step') +
  # ggsci::scale_color_jco(name = NULL) +
  scale_color_viridis_d(name = NULL) +
  theme_pubr(18) +
  theme(legend.position = 'right')

################################################################################

p1 / p2 / wrap_elements(full = p3)

ggsave('2-qc.jpeg', width = 14, height = 40)


################################################################################

filtering <- list(
  'All libraries' = qc.dat,
  'Strand-specific' = qc.dat |>
    filter(library_types != 'U'),
  '1M mRNA' = qc.dat |>
    filter(library_types != 'U') |>
    filter(non_rRNA >= 1e6),
  '40pct gene' = qc.dat |>
    filter(library_types != 'U') |>
    filter(non_rRNA >= 1e6) |>
    filter(num_mapped / non_rRNA * 100 >= 40)
)

filtering$`Min 5 libs` <-
  filtering$`40pct gene` |>
  count(genome) |>
  filter(n >= 5) |>
  left_join(filtering$`40pct gene`, 'genome')


filtering |>
  map(~ sprintf('libs: %g, gnomes: %g', nrow(.x), .x$genome |> unique() |> length()))

# $`All libraries`
# [1] "libs: 949, gnomes: 15"
# $`Strand-specific`
# [1] "libs: 791, gnomes: 14"
# $`1M mRNA`
# [1] "libs: 580, gnomes: 14"
# $`40pct gene`
# [1] "libs: 461, gnomes: 11"
# $`Min 5 libs`
# [1] "libs: 458, gnomes: 9"

browser.todo <- filtering$`Min 5 libs`
write_tsv(browser.todo, '2-todo-browser.tsv')


################################################################################

todo <- read_tsv('1-todo.tsv')

motif.tax <-
  '~/OCyRS-pipeline/data/K_motif-tax.tsv' |>
  read_tsv()

################################################################################

browser.todo |>
  select(genome, sample) |>
  unique() |>
  count(genome, name = 'RNA-seq libs') |>
  left_join(
    todo |>
      semi_join(browser.todo, c('dir' = 'genome')) |>
      left_join(motif.tax, 'tax.bio') |>
      select(dir, motif) |>
      unique() |>
      count(genome = dir, name = 'motifs'),
    'genome'
  ) |>
  arrange(desc(motifs))


################################################################################


mg <-
  motif.tax |>
  select(tax.bio, motif) |>
  unique() |>
  mutate(in.genome = TRUE) |>
  left_join(
    todo |>
      select(tax.bio) |>
      unique() |>
      left_join(motif.tax, 'tax.bio') |>
      select(tax.bio, motif) |>
      unique() |>
      mutate(has.rnaseq = TRUE),
    c('motif', 'tax.bio')
  ) |>
  left_join(
    todo |>
      semi_join(browser.todo, c('dir' = 'genome')) |>
      select(tax.bio) |>
      unique() |>
      left_join(motif.tax, 'tax.bio') |>
      select(tax.bio, motif) |>
      unique() |>
      mutate(has.rnaseq.filtered = TRUE),
    c('motif', 'tax.bio')
  ) |>
  complete(motif, tax.bio) |>
  mutate_if(is.logical, replace_na, FALSE) |>
  transmute(
    tax.bio, motif,
    x = case_when(
      has.rnaseq.filtered ~ '+ RNA-seq',
      has.rnaseq ~ '+ RNA-seq, omitted',
      in.genome ~ 'in genome',
      TRUE ~ 'not present'
    ) |>
      fct_relevel('not present', 'in genome', '+ RNA-seq, omitted', '+ RNA-seq')
  )

f.tb <-
  mg |>
  count(tax.bio, x) |>
  spread(x, n, fill = 0) |>
  arrange(desc(`+ RNA-seq`), desc(`+ RNA-seq, omitted`), desc(`in genome`)) |>
  pull(tax.bio) |>
  fct_inorder() |>
  levels()
f.mot <-
  mg |>
  count(motif, x) |>
  spread(x, n, fill = 0) |>
  arrange(desc(`+ RNA-seq`), desc(`+ RNA-seq, omitted`), desc(`in genome`)) |>
  pull(motif) |>
  fct_inorder() |>
  levels()

# sort of heatmap
x1 <-
  mg |>
  mutate_at('tax.bio', fct_relevel, f.tb) |>
  mutate_at('motif', fct_relevel, f.mot) |>
  ggplot(aes(tax.bio, motif, fill = x)) +
  geom_tile() +
  scale_fill_manual(
    name = NULL,
    values = c(
    '#F0F0F0',
    'blue',
    'orange',
    'red'
    )
  ) +
  ylab('423 motifs') +
  xlab('202 genomes (taxid + bioproject)') +
  theme_pubr(18) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )


# the combined number
x2 <-
  mg |>
  count(motif, x) |>
  mutate_at('motif', fct_relevel, f.mot) |>
  mutate_at('x', fct_relevel, 'not present', 'in genome', '+ RNA-seq, omitted', '+ RNA-seq') |>
  ggplot(aes(motif, n, fill = x)) +
  geom_bar(stat = 'identity', position = 'stack', width = 1) +
  scale_x_discrete(breaks = NULL, position = 'right') +
  scale_fill_manual(
    name = NULL,
    values = c(
    '#F0F0F0',
    'blue',
    'orange',
    'red'
    )
  ) +
  ylab('No. genomes') +
  xlab(NULL) +
  # ylim(202, 0) +
  coord_flip() +
  theme_pubr(18)

mg |>
  filter(x == '+ RNA-seq') |>
  count(motif, name = 'n.rnaseq') |>
  count(n.rnaseq) -> foo
foo |> pull(n) |> sum() -> i
x3 <-
  foo |>
  ggplot(aes(n.rnaseq, n)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = n), size = 6, nudge_y = 5) +
  xlab('No. genomes with RNA-seq') +
  ylab('No. motifs') +
  scale_x_continuous(breaks = 1:16) +
  theme_pubr(18) +
  ggtitle(sprintf('%g of 423 motifs covered by RNA-seq', i))

p1 <-
  guide_area() / (
  (x1 + theme(legend.position = 'hide')) + x2) +
  plot_layout(heights = c(1, 20)) +
  plot_layout(guides = 'collect')

((wrap_elements(full = p1) + theme_pubr(18)) / x3) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(heights = c(2, 1))

ggsave('2-motif-datasets.jpeg', width = 8, height = 12)

################################################################################