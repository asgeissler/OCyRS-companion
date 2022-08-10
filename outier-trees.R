# Create pretty phylogenetic trees for the top outliers
library(tidyverse)
library(ape)
library(ggtree)

#palette using grey
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# The top 5 outliers
xs <-  c('K07478', 'K00560', 'K06894', 'K05592', 'K01776')
names(xs) <- sprintf(
  '%s (%s)',
  c(
    'putative ATPase',
    'thymidylate synthase',
    'alpha-2-macroglobulin',
    'ATP-dependent RNA helicase DeaD',
    'galutamate racemase'
  ),
  xs
)


# Taxa data
base <- '/home/projects/rth/co2capture/subprojects/OCyRS/OCyRS-pipeline/data/'
file.path(base, 'A_representatives', 'taxonomy.tsv') %>%
  read_tsv() %>%
  mutate_at(
    c('order', 'family', 'genus', 'species'),
    str_remove,
    '^[0-9]* '
  ) %>%
  mutate_at('order', replace_na, 'unknown') -> taxa

taxa %>%
  select(order) %>%
  unique %>%
  arrange(order) %>%
  # mutate(cl = ggsci::pal_jco()(n())) %>%
  # mutate(cl = ggsci::pal_ucscgb()(n())) %>%
  mutate(cl = RColorBrewer::brewer.pal(n(), 'Paired')) %>%
  with(set_names(cl, order)) -> taxa.cl
  

# Load trees and plot colored per order
my.plot <- function(i) {
  i <- file.path(base, 'C_shrink', i, 'output.tree')
  t <- read.tree(i)
  
  t.taxa <- tibble(
    tax.bio.gca = t$tip.label
  ) %>%
    mutate(tax.bio = str_remove(tax.bio.gca, '\\.[^.]*$')) %>%
    left_join(taxa, 'tax.bio') %>%
    mutate(
      gca =  str_remove(tax.bio.gca, '^.*\\.'),
      nice = sprintf('%s (%s)', species, gca)
    )
  
  # t$tip.label <- t.taxa$nice
  t.taxa %>%
    group_by(order) %>%
    do(i = list(.$tax.bio.gca)) %>%
    # do(i = list(.$nice)) %>%
    ungroup %>%
    with(set_names(i, order)) %>%
    map(1) -> grp
  
  # p <- ggtree(t, layout = 'rectangular', branch.length = 'none')
  p <- ggtree(t, layout = 'circular', branch.length = 'none')
  groupOTU(p, grp, 'Order') +
    # geom_tiplab(as_ylab = TRUE) +
    aes(color = Order) +
    scale_color_manual(values = taxa.cl)
}

xs %>%
  map(my.plot) %>%
  invoke(
    .f = cowplot::plot_grid,
    labels = names(xs),
    hjust = 0,
    ncol = 2,
    label_size = 12
  )

ggsave('outier-trees.jpeg',
       width = 10, height = 12, dpi = 400, scale = 0.8)
