# Adapted from E_filter.R

library(tidyverse)
library(furrr)
library(ape)
library(Biostrings)

in.path <- '~/remote-server/projects/rth/co2capture/subprojects/OCyRS/OCyRS-pipeline/data/D_search-seqs-aln/*.fna.gz'
out.dir <- '~/remote-server/projects/rth/co2capture/subprojects/OCyRS/OCyRS-pipeline/data/E_search-filtered/'
out.gappy <- '~/remote-server/projects/rth/co2capture/subprojects/OCyRS/OCyRS-pipeline/data/E_search-gaps.png'

cpus <- 4
plan(multisession, workers = cpus) 

THRESHOLD <- .9

################################################################################
# Load alignments

alns <- in.path %>%
  Sys.glob() %>%
  set_names(. %>% basename %>% str_remove('.fna.gz$')) %>%
  future_map(readDNAMultipleAlignment)

################################################################################
# steps in which to count proportion of gaps
xs <- seq(0, 1, .01)

alns %>%
  future_map2(names(.), function(aln, i){
    # i <- 'K00010_downstream'
    # aln <- alns[[i]]
    # the proportion of gaps in the columns
    gap.prop <- consensusMatrix(aln, as.prob = TRUE)['-', ]
    # proportion of columns with prop. of gaps above a cutoff
    col.prop <- cut(gap.prop, xs, include.lowest = TRUE) %>%
      table %>%
      cumsum %>%
      `/`(ncol(aln))
    tibble(
      region = i,
      gappy.threshold = xs[-1] * 100,
      col.prop = col.prop * 100
    )
  }) %>%
  bind_rows() -> gappy

gappy %>%
  ggplot(aes(gappy.threshold, col.prop, group = region)) +
  geom_line(alpha = 0.2, color = 'black') +
  geom_smooth(aes(group = 1), color = 'blue') +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  geom_vline(xintercept = 90, color = 'red') +
  ylab('Proportion of alignment columns %') +
  xlab('Threshold % gap positions per column') +
  theme_bw(18)

ggsave('supplement-figure_search-gaps.png', width = 8, height = 6)
