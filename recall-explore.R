# Build potential hypothesis on why some Rfam families have not been recalled

# - quality of Rfam hit (partial hits?)
# - No hits vs overlaps

# re-uses code from J_recall.R in pipeline

library(tidyverse)

library(plyranges)
library(conflicted)

conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("n", "dplyr")
conflict_prefer("first", "dplyr")


################################################################################

in.cats <- '../OCyRS-pipeline/data/J_FDR-categories.tsv'
in.ref <- '../OCyRS-pipeline/data/J_novel/references_inside-of_intergenic_regions.tsv.gz'


# Colorblind-friendly palettes of the Color Universal Design
# https://riptutorial.com/r/example/28354/colorblind-friendly-palettes
cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7")

################################################################################

cats <- read_tsv(in.cats)
ref <- read_tsv(in.ref)


ranges <- as_granges(ref)

################################################################################
# Compute overlaps, but only per group to prevent erroneous multiple counts

ranges %>%
  filter(type == 'Candidate motif') %>%
  reduce_ranges_directed() %>%
  mutate(reduced.row = 1:plyranges::n()) -> ranges2.red.motifs

ranges %>%
  filter(type != 'Candidate motif') %>%
  group_by(type, rfam) %>%
  reduce_ranges_directed() %>%
  ungroup %>%
  mutate(reduced.row = 1:plyranges::n()) -> ranges2.red.ref

join_overlap_intersect(
  ranges2.red.motifs,
  ranges2.red.ref %>%
    mutate(strand2 = strand)
) %>%
  mutate(
    overlap = IRanges::width(.),
    orientation = ifelse(
      strand == strand2,
      'sense',
      'anti-sense'
    )
  ) -> overlaps


################################################################################
# stats on overlaps relative to reduced positions

over.view <-
  ranges2.red.ref |>
  as_tibble() |>
  count(type, rfam, name = 'nredhits') |>
  left_join(
    overlaps |>
      as_tibble() |>
      select(type, rfam, reduced.row.y) |>
      unique() |>
      count(type, rfam, name = 'overlapred'),
    c('type', 'rfam')
  ) |>
  mutate_at('overlapred', replace_na, 0)


################################################################################

ref.overview <-
  ref |>
  filter(type != 'Candidate motif')  |>
  count(type, name, rfam, region) |>
  group_by(type, name, rfam) |>
  summarize(
    nregions = n(),
    nhits = sum(n),
    avghits = mean(n)
  ) |>
  ungroup() |>
  left_join(
    ref |>
      mutate(tax.bio = str_remove(seqnames, '\\.[^.]*$')) |>
      select(type, name, rfam, tax.bio) |>
      unique() |>
      count(type, name, rfam, name = 'species'),
    c('type', 'name', 'rfam')
  ) |>
  left_join(over.view, c('type', 'rfam'))

################################################################################

p1 <-
  ref.overview |>
  filter(type == 'Rfam') |>
  mutate(r = overlapred / nredhits) |>
  ggplot(aes(nredhits, species)) +
  geom_point() +
  annotation_logticks(side = 'b') +
  scale_x_log10(labels = scales::comma) +
  xlab('CMsearch hits') +
  theme_bw(18)

p2 <-
  ref.overview |>
  filter(type == 'Rfam') |>
  mutate(r = overlapred / nredhits) |>
  ggplot(aes(nredhits, overlapred)) +
  geom_point() +
  xlab('CMsearch hits') +
  ylab('CMsearch hits with motif overlap\n(sense and/or anti-sense)') +
  scale_x_log10(labels = scales::comma) +
  annotation_logticks(side = 'b') +
  theme_bw(18)

p3 <-
  ref.overview |>
  filter(type == 'Rfam') |>
  group_by(nredhits) |>
  summarize(
    fams = n(),
    recall = sum(overlapred > 0)
  ) |>
  ungroup() |>
  arrange(desc(nredhits)) |>
  mutate(
    cumfam = cumsum(fams),
    cumrecall = cumsum(recall)
  ) |>
  ggplot(aes(nredhits)) +
  geom_line(aes(y = cumfam), color = 'blue') +
  geom_line(aes(y = cumrecall / cumfam * 55), color = 'red') +
  scale_y_continuous(
    # Features of the first axis
    name = "No. Rfam families",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~ . / 55, name=" Rate of family recall")
  ) +
  xlab('Minimal no. CMsearch hits') +
  scale_x_log10(labels = scales::comma) +
  annotation_logticks(side = 'b') +
  theme_bw(18) +
  theme(
    axis.title.y = element_text(color = 'blue'),
    axis.title.y.right = element_text(color = 'red')
  )


################################################################################

library(patchwork)

(p1 + p2) / p3 +
  plot_annotation(tag_levels = 'A')

ggsave('recall-explore.jpeg', width = 12, height = 10)


################################################################################
# species/hits plot, but only for search region with most hits per family

# ref |>
#   dplyr::filter(type == 'Rfam') |>
#   count(name, region) |>
#   group_by(name) |>
#   slice_max(n) |>
#   ungroup() |>
#   select(- n) |>
#   left_join(ref, c('name', 'region')) -> ref2
# ref <-
#   ref |>
#   filter(type != 'Rfam') |>
#   bind_rows(ref2)

# re-run above; did not really change anything
# -> not enough hits for CMfinder to detect
