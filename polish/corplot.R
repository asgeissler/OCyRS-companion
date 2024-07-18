# Adapted from J_categorize
# Inspect distribution of scores
# Further classify by variation and power

library(tidyverse)
library(corrplot)

in.fdr <- '~/remote-server/projects/rth/co2capture/subprojects/OCyRS/OCyRS-pipeline/data/I_fdr.tsv'
in.scores <- '~/remote-server/projects/rth/co2capture/subprojects/OCyRS/OCyRS-pipeline/data/H_scores.tsv'
in.cmstat <- '~/remote-server/projects/rth/co2capture/subprojects/OCyRS/OCyRS-pipeline/data/I_cmstat.tsv'

fdr <- read_tsv(in.fdr)
scores <- read_tsv(in.scores)
cmstat <- read_tsv(in.cmstat)


###############################################################################
fdr <- read_tsv(in.fdr)

  
fdr %>%
  filter(RNAphylo.fdr <= 10) %>%
  select(motif) %>%
  mutate(
    dir = 'D_search-seqs',
    cat = 'FDR ≤ 10%'
  ) -> cats1

###############################################################################

score.dat <- scores %>%
  left_join(
    cmstat %>%
      select(motif, bifurcations, rel_entropy_cm, rel_entropy_hmm) %>%
      mutate(dir = ifelse(
        str_starts(motif, 'RF'),
        'G_rfam-bacteria-seeds',
        'D_search-seqs'
      )),
    c('motif', 'dir')
  ) %>%
  transmute(
    motif, dir,
    # CMfinder pipeline stats
    'Length, log10' = log10(alen + 1),
    'Sequences, log10' = log10(nseq + 1),
    'RNAphylo, log10' = log10(RNAphylo + 1),
    'hmmpair, log10' = log10(hmmpair + 1),
    # Ratio
    # 'log10(RNAphylo / Length)' = log10(RNAphylo / alen),
    # R-Scape
    'Average SI %' = avgid,
    'Paired positions %' = 2 * nbpairs / alen * 100,
    'Alignment power %' = expected / nbpairs * 100,
    'Covarying bps %' = observed / nbpairs * 100,
    # CMstat
    # consensus_residues_len, expected_max_hit_len,
    # Bifurcations = bifurcations,
    # 'Rel. entropy CM' = rel_entropy_cm,
    # 'Rel. entropy hmm' =  rel_entropy_hmm
  )

###############################################################################



pdf('fig_FDR-score-correlations.pdf')
cats1 %>%
  left_join(score.dat, c('dir', 'motif')) %>%
  select(- c(motif, dir, cat))  %>%
  # GGally::ggpairs()
  cor %>%
  corrplot(
    method = 'square',
    order = 'AOE',
    type = 'lower',
    diag = FALSE,
    insig='blank',
    addCoef.col ='black',
    number.cex = 0.8
  )
dev.off()

