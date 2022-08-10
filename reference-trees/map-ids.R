library(tidyverse)
library(ape)

library(rentrez)

# Purpose: Map leaves in nodes to taxonid


################################################################################
################################################################################
# Step 1: Load individual trees and rename tips to genbank/refseq values
#         Only if possible, of course.

#-------------------------------------------------------------------------------
# Cornet
orper <- read.tree('2022_09_06-ref-trees/Vertical-ORPER.tre')

orper$tip.label %>%
  keep(str_detect, '_GC[AF]_') -> with.gcf

with.gcf %>%
  str_remove('^.*(?=GC[FA]_)') %>%
  set_names(with.gcf) -> gcf

orper2 <- keep.tip(orper, with.gcf)
orper2$tip.label <-  gcf[orper2$tip.label]

#-------------------------------------------------------------------------------
# Chen, on marker
chen.marker <- read.tree('chen-amphora-low.txt')

chen.marker$tip.label %>%
  keep(str_detect, '_GC[AF]_') -> with.gcf

with.gcf %>%
  str_remove('^.*(?=GC[FA]_)') %>%
  set_names(with.gcf) -> gcf

chen.marker2 <- keep.tip(chen.marker, with.gcf)
chen.marker2$tip.label <-  gcf[chen.marker2$tip.label]

#-------------------------------------------------------------------------------
# Chen, on genome, generally better
chen.busco <- read.tree('chen-busco-high.txt')

chen.busco$tip.label %>%
  keep(str_detect, '_GC[AF]_') -> with.gcf

with.gcf %>%
  str_remove('^.*(?=GC[FA]_)') %>%
  set_names(with.gcf) -> gcf

chen.busco2 <- keep.tip(chen.busco, with.gcf)
chen.busco2$tip.label <-  gcf[chen.busco2$tip.label]
#-------------------------------------------------------------------------------
# Moore

moore <- read.nexus('moore-Figure1.tre')
# only in scientific name format

################################################################################
# For trees with GCA/GCF lookup taxid
gca.trees <- list(
  'Cornet' = orper2,
  'Chen_multigene' = chen.marker2,
  'Chen_genomic' = chen.busco2
)

# The assemblies used
gca.trees %>%
  map('tip.label') %>%
  unlist() %>%
  unique -> gcf

# chunk into smaller pieces for query uid of assembly
split(
  gcf,
  ceiling(seq_along(gcf) / 10)
) %>%
  map(str_c, collapse = ' OR ') %>%
  map(~ entrez_search('assembly', term = .x)) %>%
  map('ids') %>%
  unlist -> gcf.uids

# Speed up linking with web history, though even that is limited
# chunk again in 100er block
split(
  gcf.uids,
  ceiling(seq_along(gcf.uids) / 100)
) %>%
  map(function(gcf.i) {
    whist <- entrez_post('assembly', id = gcf.i)
    
    # establish link assembly <-> taxonomy
    lnk <- entrez_link(dbfrom = 'assembly', db = 'taxonomy', web_history = whist)
    whist.tax <- entrez_post('taxonomy', id = lnk$links$assembly_taxonomy)
    
    # match GCF, taxonid, and species name
    sum.asm <- entrez_summary('assembly', web_history = whist)
    sum.tax <- entrez_summary('taxonomy', web_history = whist.tax)
    
    # process summaries to lookup table to rename tree
    sum.tax %>%
      map(as_tibble) %>%
      map(select, taxid, scientificname) %>%
      bind_rows() -> tax
    sum.asm %>%
      map(extract_from_esummary, c('synonym', 'biosampleaccn', 'taxid')) %>%
      map(as_tibble) %>%
      bind_rows() %>%
      unnest(synonym) %>%
      filter(synonym != 'identical') %>%
      mutate_at('taxid', as.integer) %>%
      left_join(tax, 'taxid') %>%
      mutate_at('scientificname', str_replace_all, '[^A-Za-z0-9]', '-') %>%
      mutate(taxbioname = paste(taxid, biosampleaccn, scientificname, sep = '.'))
  }) %>%
  bind_rows() -> look

look %>%
  filter(synonym != '') %>%
  with(set_names(taxbioname, synonym)) -> look2

# Rename tips in the gca.trees
gca.trees %>%
  map(function(i) {
    i$tip.label <- look2[i$tip.label]
    i
  }) -> nice.trees

################################################################################

moore$tip.label %>%
  str_trim() -> qs

qs %>%
  split(ceiling(seq_along(qs) / 10)) %>%
  map(str_c, collapse = ' OR ') %>%
  map(function(q) {
    res <- entrez_search('taxonomy', q)
    
    whist <- entrez_post('taxonomy', id = res$ids)
    entrez_fetch('taxonomy', web_history = whist,
                 rettype = 'xml', retmax = 200)
  }) -> res

res %>%
  map(function(i) {
    res.x <- xml2::read_xml(i)
    
    xml2::xml_children(res.x) %>%
      map(function(i) {
        # i <- xml2::xml_children(res.x)[[1]]
        
        list(
          taxid =  xml2::xml_find_all(i, 'TaxId') %>%
            xml2::xml_text(),
          name =  xml2::xml_find_all(i, 'ScientificName') %>%
            xml2::xml_text(),
          syn = xml2::xml_find_all(i, 'OtherNames/Synonym') %>%
            map(xml2::xml_text) %>%
            unlist,
          disp = xml2::xml_find_all(i, 'OtherNames/Name/DispName') %>%
            map(xml2::xml_text) %>%
            unlist
        )
      })
  }) %>%
  reduce(c) -> res2

res2 %>%
  map(function(x) {
    x$ns <- c(x$name, x$syn, x$disp)
    x$name <- NULL
    x$syn <- NULL
    x$disp <-NULL
    x
  }) %>%
  map(as_tibble) %>%
  bind_rows() %>%
  mutate(
    nice = paste(
      taxid,
      str_replace_all(ns, '[^A-Za-z0-9]+', '-'),
      sep = '.'
    ),
    ns = str_replace_all(ns, '[ .-]+', '_')
  ) %>%
  unique -> tax
  
  
mask <- qs %in% tax$ns
table(mask)
# 74.1% -> great

with(tax, set_names(nice, ns)) -> look

moore2 <- keep.tip(moore, moore$tip.label[moore$tip.label %in% names(look)])
moore2$tip.label <- look[moore2$tip.label]

nice.trees$Moore <- moore2

################################################################################
# Save out nicer trees

map2(
  sprintf('nice-trees/%s.tree', names(nice.trees)),
  nice.trees,
  ~ write.tree(.y, file = .x)
)
