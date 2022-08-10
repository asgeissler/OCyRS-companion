library(tidyverse)
library(knitr)
library(magrittr)
library(cowplot)
library(gridExtra)

library(gt)

library(ape)
library(ggtree)
library(tidytree)
library(patchwork)

nw <- '((A:1,B:0.7):0.3,((C:1,D:0.8):0.4,(E:0.2, F:0.3):0.4):0.1):0;'
nw2 <- '((A:1,B:0.7):0.5,((E:1,D:0.8):0.4,C:0.6):0.1):0;'

t1 <- read.tree(text = nw)
t2 <- read.tree(text = nw2)

t1 %>%
  as.treedata() %>%
  mutate(hl = ifelse(label == 'F', 'red', 'blue')) %>%
  ggtree(ladderize = FALSE) +
  geom_tiplab(aes(color = I(hl)), show.legend = FALSE) +
  geom_label(aes(x = branch, label = branch.length)) +
  theme_tree2() +
  ggtitle('Example tree t1') -> p1
  
t1 %>%
  drop.tip('F') -> t1b
t1b %>%
  ggtree(ladderize = FALSE) +
  geom_tiplab(color = 'blue') +
  geom_label(aes(x = branch, label = branch.length)) +
  theme_tree2() +
  ggtitle("t1' = drop.tip(t1, 'F')") -> p2

t2 %>%
  as.treedata() %>%
  mutate(hl = ifelse(label %in% c('E', 'C'), 'red', 'blue')) %>%
  ggtree(ladderize = FALSE) +
  geom_tiplab(aes(color = I(hl)), show.legend = FALSE) +
  geom_label(aes(x = branch, label = branch.length)) +
  theme_tree2() +
  ggtitle('Example tree t2') -> p3

tribble(
  ~ split, ~ f_t1, ~ f_t2,
  'AB,CDE', 0.4, 0.6,
  'CD,ABE', 0.4, 0,
  'DE,ABC', 0, 0.4,
) |>
  mutate(
    diff = abs(f_t1 - f_t2),
    squared.diff = diff ** 2
  ) |>
  mutate_at('diff', as.character) %>%
  add_row(diff = 'Sum:', squared.diff = sum(.$squared.diff)) %>%
  add_row(diff = 'Sqrt:', squared.diff = sqrt(last(.$squared.diff)) |> round(2)) |>
  mutate_all(as.character) |>
  mutate_all(replace_na, '') |>
  tableGrob(rows = NULL, theme = ttheme_minimal(base_size = 18)) -> p4

(p1 | p2 | p3) / p4

ggsave('supplement-fig_topo.jpg',
       scale = .7,
       width = 12, height = 8)


assertthat::are_equal(
  dist.topo(t1b, t2, method = 'score') %>%
    as.matrix %>% `[`(1,2),
  sqrt(.36),
  tol = .Machine$double.eps
)
