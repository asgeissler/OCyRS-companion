# sanitize intuition, any CRS predicitons overlapping Brenes predicted sRNA?

library(tidyverse)
library(plyranges)

################################################################################

srna <-
  '2024-04-08-Brenes-Supplement.XLSX' |>
  xlsx::read.xlsx(sheetName = 'Total', startRow = 9 ) |>
  transmute(
    seqnames = '103690.SAMD00061094.BA000019',
    start = pmin(TSS, End),
    end = pmin(TSS, End),
    strand = ifelse(Strand == 'fwd', '+', '-'),
    srna = 1:n()
  ) |>
  as_granges() 

################################################################################

mpos <-
  '/home/projects/rth/co2capture/subprojects/OCyRS/OCyRS-pipeline/data/J_motif-aln-seq-pos.tsv' |>
  read_tsv() |>
  filter(chr == '103690.SAMD00061094.BA000019') |>
  transmute(
    seqnames = chr,
    start,
    end,
    strand,
    motif,
    id = paste(motif, start)
  ) |>
  as_granges() 

################################################################################

novel.crs <-
  '/home/projects/rth/co2capture/subprojects/OCyRS/OCyRS-pipeline/data/J_novel/potentially-novel-motifs.tsv' |>
  read_tsv()

novel.red <-
  '/home/projects/rth/co2capture/subprojects/OCyRS/OCyRS-pipeline/data/L_redundant.tsv' |>
  read_tsv()

novel.ta <-
  '../Public-RNAseq/5-expression-ratios.tsv' |>
  read_tsv()


################################################################################

helper <- function(x, mpos) {
  foo <- mpos |> join_overlap_intersect(srna)
  
  c(
    motifs =  mpos |>
      as_tibble() |>
      pull(motif) |>
      unique() |>
      length(),
    motif.pos = mpos |> length(),
    motif.overlap =  foo |>
      as_tibble() |>
      pull(motif) |>
      unique() |>
      length(),
    srna.overlap =  foo |>
      as_tibble() |>
      pull(srna) |>
      unique() |>
      length(),
    x = x
  )
}
  

################################################################################

list(
  'FDR max 10% CMfinder predictions' = mpos$motif |> unique(),
  'Novel CRS' = novel.crs |>
    pull(motif) |>
    intersect(mpos$motif),
  'Novel CRS, non-redundant' = novel.crs |>
    anti_join(novel.red, c('motif' = 'redundant')) |>
    pull(motif) |>
    intersect(mpos$motif),
  'Novel CRS, non-redundant, transcriptionally active in one genome' = novel.crs |>
    anti_join(novel.red, c('motif' = 'redundant')) |>
    semi_join(
      novel.ta |>
        filter(expressed.libs >= 3) |>
        mutate(motif = str_remove(gene, ';pos.*$')),
      'motif'
    ) |>
    pull(motif) |>
    intersect(mpos$motif),
  'Novel CRS, non-redundant, transcriptionally active in PCC 7120' = novel.crs |>
    anti_join(novel.red, c('motif' = 'redundant')) |>
    semi_join(
      novel.ta |>
        filter(genome == 'Nostoc.sp.PCC.7120.FACHB.418_txid103690_SAMD00061094_Nostocales') |>
        filter(expressed.libs >= 3) |>
        mutate(motif = str_remove(gene, ';pos.*$')),
      'motif'
    ) |>
    pull(motif) |>
    intersect(mpos$motif)
) |>
  map(unique) |>
  map(~ filter(mpos, motif %in% .x)) %>%
  map2(names(.), ~ helper(.y, .x)) |>
  bind_rows() -> res

res |>
  select(x, everything())

# x                                                                motifs motif.pos motif.overlap srna.overlap
# 1 FDR max 10% CMfinder predictions                                 139    201       10            7           
# 2 Novel CRS                                                        36     49        0             0           
# 3 Novel CRS, non-redundant                                         35     48        0             0           
# 4 Novel CRS, non-redundant, transcriptionally active in one genome 26     28        0             0           
# 5 Novel CRS, non-redundant, transcriptionally active in PCC 7120   20     20        0             0     