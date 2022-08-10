# prepare RNA-Schlange pipelines for downloading and pre-processing the
# selected RNA-seq datasets (≥ 5 libraries)

library(tidyverse)

################################################################################

srp <-
  '0-public-rnaseq-projects.tsv' |>
  read_tsv()

over <-
  '0-public-stat.tsv' |>
  read_tsv()

################################################################################
# Selection: ≥5 libraries

dat.selection <-
  over |>
  filter(no.libs >= 5) |>
  select(tax, order, species) |>
  unique() |>
  mutate_if(is.character, str_replace_all, '[^A-Za-z0-9]', '.') |>
  left_join(
    srp |>
      filter(tax == tax.search),
    'tax'
  )

################################################################################
# Inconvenience: 1 taxid might match multiple bioprojects
# (tax+biorpoject is the unique identifier of genome in proGenomes)
# -> match and check if running multiple times if needed

progenomes <-
  '~/OCyRS-pipeline/data/A_representatives/taxonomy.tsv' |>
  read_tsv() |>
  select(tax.bio) |>
  separate(tax.bio, c('tax', 'bio'), sep = '\\.', remove = FALSE)

todo <-
  dat.selection |>
  select(tax) |>
  unique() |>
  mutate_at('tax', as.character) |>
  left_join(progenomes, 'tax') |>
  # count(tax) |> count(n)
  # n    nn
  # 1     1    17
  # ok, multiple runs not needed, but still keep bioproject in dir name for
  # more transparency
  mutate_at('tax', as.integer) |>
  left_join(dat.selection, 'tax') |>
  transmute(
    SRA, SRR, SRP,
    tax.bio,
    dir = sprintf(
      '%s_txid%s_%s_%s',
       species, tax, bio, order
    ) |>
      str_replace_all('\\.\\.+', '.')
  ) |>
  separate_rows(SRR, sep = ';')

write_tsv(todo, '1-todo.tsv')


################################################################################
# Create one RNA-Schlange project per folder and export the SRA list for 
# downloading and processing

dir.create('1-RNA-Schlange')

todo |>
  split(todo$dir) |>
  map(function(xs) {
    # xs <- filter(todo, dir == 'Nostoc.sp.PCC.7120.FACHB.418_txid103690_SAMD00061094_Nostocales')
    x <- first(xs$dir)
    taxbio <- first(xs$tax.bio)
    # Create workflow for dir
    cmd <- sprintf(
      'git clone git@github.com:asgeissler/RNA-Schlange.git 1-RNA-Schlange/%s',
      x
    )
    system(cmd)
    
    # export SRA info
    xs |>
      transmute(
        run = SRR,
        condition = paste(
          SRP, SRA,
          sep = '-'
        )
      ) |>
      write_csv(sprintf(
        '1-RNA-Schlange/%s/sra-SE.csv',
        x
      ))
    
    # copy genome over from main project
    file.copy(
      '~/OCyRS-pipeline/data/A_representatives/%s/genome.fna.gz' |>
        sprintf(taxbio),
      '1-RNA-Schlange/%s/genome.fna.gz' |>
        sprintf(x)
    )
    
    
    # make an empty annotation file but use the sequences instead for
    # the salmon quantification
    file.create(sprintf('1-RNA-Schlange/%s/genome.gff.gz', x))
    dir.create(sprintf('1-RNA-Schlange/%s/analysis', x))
    file.copy(
      '~/OCyRS-pipeline/data/A_representatives/%s/genes.fasta.gz' |>
        sprintf(taxbio),
      '1-RNA-Schlange/%s/analysis/30_genes.fna.gz' |>
        sprintf(x)
    )
  })

# now run:
# bash 1-donwload.sh

# TODO full run script
# Then run the individual jobs on the cluster
# bash 1-full-slurm.sh

################################################################################
# manually confirm that all the downloads were successful and have the right
# SRR download files (not empty)

todo <- read_tsv('1-todo.tsv') |>
  # misclassified PE
  filter(SRP != 'SRP371565') |>
  # not normal RNAseq
  filter(SRR != 'SRR2079408') |>
  filter(SRR != 'SRR24734921') |>
  filter(SRR != 'SRR24734922')
foo <-
  todo |>
  pull(dir) |>
  unique() |>
  map(function(i) {
    # i <- 'Prochlorococcus.marinus.str.MIT.9312_txid74546_SAMN02598321_NA'
    # i <- 'Synechocystis.sp.PCC.6803_txid1148_SAMD00061113_Synechococcales'
    expected.srr <-
      file.path('1-RNA-Schlange', i, 'sra-SE.csv') |>
      read_csv() |>
       pull(run) |>
       sort()
    
    downloaded.srr <-
      file.path('1-RNA-Schlange', i, 'data/*fastq.gz') |> 
      Sys.glob() |>
      basename() |>
      str_remove('.fastq.gz') |>
      sort()
    
    
    tibble(
      dir = i,
      error.donwload.log = file.path('1-RNA-Schlange', i, 'log-download.err.txt') |>
        read_lines() |>
        str_detect('[Ee][Rr][Rr][Oo][Rr]') |>
        any(),
      all.downloaded =  identical(downloaded.srr, expected.srr),
      full.run = file.path('1-RNA-Schlange', i, 'analysis/53_main_multiqc/multiqc_data/multiqc_data.json') |>
        file.exists()
    )
  }) |>
  bind_rows()

foo |>
  filter(all.downloaded, !full.run) |>
  pull(dir) |>
  str_c(collapse = ',')

# bash run_slurm_singularity.sh  --set-resources ribo_removal:mem_mb=32GB  ribo_removal:time='10:00:00' --cores all all

################################################################################

# datasets downloaded
foo |>
  count(all.downloaded) |>
  spread(all.downloaded, n) |>
  (\(x) if (! 'FALSE' %in% colnames(x)) mutate(x, 'FALSE' = 0) else  x)() |>
  mutate(x = `TRUE`, y = `TRUE` + `FALSE`) |>
  with(sprintf('%g of %g genomes (%.2f%%)', x, y, x / y * 100))

# % libs downloaded
foo |>
  left_join(todo) |>
  count(all.downloaded) |>
  spread(all.downloaded, n) |>
  (\(x) if (! 'FALSE' %in% colnames(x)) mutate(x, 'FALSE' = 0) else  x)() |>
  mutate(x = `TRUE`, y = `TRUE` + `FALSE`) |>
  with(sprintf('%g of %g libraries (%.2f%%)', x, y, x / y * 100))
# After removing libraries without any QC passing reads:
# "908 of 1036 libraries (87.64%)"
  

# datasets completed
foo |>
  count(full.run) |>
  spread(full.run, n) |>
  (\(x) if (! 'FALSE' %in% colnames(x)) mutate(x, 'FALSE' = 0) else  x)() |>
  mutate(x = `TRUE`, y = `TRUE` + `FALSE`) |>
  with(sprintf('%g of %g datasets (%.2f%%)', x, y, x / y * 100))

################################################################################
# libraries with exactly 0 QC passing reads

qc.dat <-
  '1-RNA-Schlange/*/analysis/51_fastqc_after/multiqc_data/multiqc_fastqc.txt' |>
  Sys.glob() %>%
  set_names(., .) |>
  map(read_tsv)

qc.dat |>
  map(nrow) |>
  unlist() |>
  sum() -> libs.total

qc.dat |>
  map(filter, `Total Sequences` > 0) |>
  map(nrow) |>
  unlist() |>
  sum() -> libs.with.reads

sprintf(
 '%g of %g libs have reads after mapping (%.1f%%)',
  libs.with.reads,
  libs.total,
  libs.with.reads / libs.total * 100
)

qc.failed <-
  qc.dat |>
  map(filter, `Total Sequences` == 0) |>
  keep(~ nrow(.x) > 0)



# samples to remove
qc.rm <-
  set_names(
    qc.failed |>
      map(pull, Sample),
    qc.failed |>
      names() |>
      str_remove('^1-RNA-Schlange/') |>
      str_remove('/.*$')
  )

# backup the original samples
qc.rm |>
  names() |>
  sprintf(fmt = '1-RNA-Schlange/%s/samples.csv') %>%
  file.copy(., paste0(., '.bak'))

rm.helper <- function(qc.rm) {
  # remove offending lines
  qc.rm |>
    names() |>
    sprintf(fmt = '1-RNA-Schlange/%s/samples.csv') %>%
    set_names(., .) |>
    map(read_csv) |>
    map(mutate, Sample = paste('batch1', sample, condition, sep = '_')) |>
    map2(qc.rm, ~ filter(.x, ! (Sample %in% .y))) |>
    map(select, - Sample) %>%
    map2(names(.), ~ write_csv(.x, .y))
  
  # move fastq from data to backup folder
  qc.rm |>
    names() |>
    sprintf(fmt = '1-RNA-Schlange/%s/data.bad.qc') |>
    map(dir.create)
  
  qc.rm |>
    map(str_remove, '^batch1_') |>
    map(str_remove, '_.*$') %>%
    map2(names(.), ~ sprintf('1-RNA-Schlange/%s/data.bad.qc/%s.fastq.gz', .y, .x)) |>
    unlist(use.names = FALSE) |>
    map(~ file.rename(str_remove(.x, '.bad.qc'), .x))
}

rm.helper(qc.rm)

################################################################################
# clear bad libraries that don't match reference

# error <- 'salmon was only able to assign 0 fragments to transcripts in the index, but the minimum number of required assigned fragments (--minAssignedFrags) was 10. This could be indicative of a mismatch between the reference and sample, or a very bad sample.'
error <- 'salmon was only able to assign 0 fragments to transcripts in the index'

salmon.qc <-
  '1-RNA-Schlange/*/logs/salmon/*.log' |>
  Sys.glob() %>%
  set_names(., .) |>
  map(read_lines) |>
  keep(~ .x |> str_detect(error) |> any())

salmon.rm <-
  tibble(i = names(salmon.qc)) |>
  mutate(
    org = i |>
      str_remove('^1-RNA-Schlange/') |>
      str_remove('/.*$'),
    x = i |>
      basename() |>
      str_remove('.log$')
  ) |>
  group_by(org) |>
  do(xs = list(.$x)) |>
  with(set_names(xs, org)) |>
  map(1)

rm.helper(salmon.rm)
