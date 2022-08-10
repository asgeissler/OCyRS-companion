# Computes the ratio of orthologous gene (OG) anchors vs
# all number of genes in the genome.
# Coverage per genome such as by median/mean and standard deviation.

library(tidyverse)

all.genes <-
  '~/OCyRS-pipeline/data/A_representatives/genes.tsv.gz' |>
  read_tsv()

ogs  <-
  '~/OCyRS-pipeline/data/B_OGs.tsv' |>
  read_tsv()


inner_join(
  ogs |>
    select(tax.bio, tax.bio.gene) |>
    unique() |>
    count(tax.bio, name = 'no.og'),
  all.genes |>
    select(tax.bio, tax.bio.gene) |>
    unique() |>
    count(tax.bio, name = 'no.genes'),
  'tax.bio'
) -> foo
  
foo |>
  mutate(prop = no.og / no.genes * 100) |>
  pull(prop) |>
#   summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 15.39   22.14   25.05   25.42   29.10   39.21 
  sd()
# [1] 4.274378
