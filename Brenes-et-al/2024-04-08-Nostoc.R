# Brenes et al did their TSS/sRNA in Nostoc sp. PCC 7120  txid103690
# ahead of their RNA prediction
# -> What are our numbers in this species?


library(tidyverse)
library(plyranges)

search <-
  '~/OCyRS-pipeline/data/J_novel/all_intergenic_regions.tsv.gz' |>
  read_tsv()


search |>
  filter(str_detect(seqnames, '^103690\\.')) |>
  head()
  # count(seqnames) -> 1134
  # count(region) |> nrow() 872
  # count(origin.gene) |> nrow() 663


################################################################################
# overlap comparison with Brenes

ref <-
  '2024-04-08-Brenes-Supplement.XLSX' |>
  xlsx::read.xlsx(sheetName = 'Total', startRow = 9 )

ref |>
  transmute(
    seqnames = '103690.SAMD00061094.BA000019',
    start = pmin(TSS, End),
    end = pmin(TSS, End),
    strand = ifelse(Strand == 'fwd', '+', '-'),
    srna = 1:n()
  ) |>
  as_granges() -> ref.range

search |>
  filter(str_detect(seqnames, '^103690\\.')) |>
  as_granges() -> search.range

join_overlap_inner(
  ref.range,
  search.range
) |>
  as_tibble() |>
  pull(srna) |>
  unique() |>
  length() -> foo


# > foo / nrow(ref) * 100
# [1] 18.96024
# > foo
# [1] 62
# > nrow(ref)
# [1] 327