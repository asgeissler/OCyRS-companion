# 1. RNA-browser git repositories with links to RNA-Schlange output
# 2. Prepare annotation for quantification
# 3. Submit jobs


library(tidyverse)


################################################################################

todo.schlange <-
  '1-todo.tsv' |>
  read_tsv()

todo.browser <-
  '2-todo-browser.tsv' |>
  read_tsv() |>
  select(dir = genome, file = sample, library_types) |>
  separate(file, c('batch', 'SRR', 'SRP', 'SRA'),
           sep = '[_-]', remove = FALSE) |>
  unique() |>
  left_join(
    todo.schlange |>
      select(tax.bio, dir) |>
      unique(),
    'dir'
  )

# motifs to exclude due to redundancy
redundant <-
  '/home/projects/rth/co2capture/subprojects/OCyRS/OCyRS-pipeline/data/L_redundant.tsv' |>
  read_tsv() |>
  select(motif = redundant) |>
  unique()


################################################################################

motif.pos <-
  '~/OCyRS-pipeline/data/K_motif-tax-pos.tsv' |>
  read_tsv() |>
  anti_join(redundant, 'motif')

genes <-
  '~/OCyRS-pipeline/data/A_representatives/genes.tsv.gz' |>
  read_tsv()

# Prepare for species of interest
mask <-
  todo.browser |>
  pull(tax.bio) |>
  unique()
# SAF columns are: 'GeneID\tChr\tStart\tEnd\tStrand'
todo.genes <-
  bind_rows(
    motif.pos |>
      filter(tax.bio %in% mask) |>
      transmute(
        tax.bio,
        GeneID = sprintf('%s;pos%s', motif, start),
        Chr = seqnames,
        Start = start,
        End = end,
        Strand = strand
      ) |>
      mutate_all(as.character),
    genes |>
      filter(tax.bio %in% mask) |>
      transmute(
        tax.bio,
        GeneID = tax.bio.gene,
        Chr = tax.bio.chr,
        Start = start,
        End = end,
        Strand = strand
      ) |>
      mutate_all(as.character)
)

################################################################################
# Create one RNA-browser project per folder and copy processed mRNA over

# pattern location to the mRNA from the pre-processing
pattern <- '1-RNA-Schlange/%s/analysis/22_non-ribosomal_RNA/%s.fastq.gz'

dir.create('3-RNA-Browser')

todo.browser |>
  split(todo.browser$dir) |>
  map(function(xs) {
    # xs <- filter(todo.browser, dir == 'Acaryochloris.marina.MBIC11017_txid329726_SAMN02604308_unknown')
    x <- first(xs$dir)
    taxbio <- first(xs$tax.bio)
    # Create workflow for dir
    cmd <- sprintf(
      'git clone git@github.com:asgeissler/RNA-Browser.git 3-RNA-Browser/%s',
      x
    )
    system(cmd)

    # copy genome over from main project
    file.copy(
      '~/OCyRS-pipeline/data/A_representatives/%s/genome.fna.gz' |>
        sprintf(taxbio),
      '3-RNA-Browser/%s/genome.fna.gz' |>
        sprintf(x)
    )

    # copy (not symbolic link due to singularity container security limitation)
    # the processed read files over
    dir.create(sprintf('3-RNA-Browser/%s/data', x))
    xs |>
      pull(file) |>
      map(function(i) {
        file.copy(
          sprintf(pattern, x, i),
          sprintf('3-RNA-Browser/%s/data/%s.fastq.gz', x, i)
        )
      })

    # create the samples file
    xs |>
      transmute(
        dataset = SRP,
        reverse = ifelse(library_types == 'SR', 'yes', 'no'),
        sample = SRA,
        condition = SRR,
        file = paste0(file, '.fastq.gz')
      ) |>
      write_csv(sprintf('3-RNA-Browser/%s/samples.csv', x))

    # make an empty annotation file
    file.create(sprintf('3-RNA-Browser/%s/genome.gff.gz', x))
    # build SAF for featureCounts
    todo.genes |>
      filter(tax.bio == taxbio) |>
      select(- tax.bio) |>
      write_tsv(sprintf('3-RNA-Browser/%s/genes.saf', x))
  })

# rm -rf 3-RNA-Browser/*/analysis/4*
# update config files to
# vi 3-RNA-Browser/*/config.yaml
# '--fracOverlap 0.5 -O'
# for the featureCount options
# double check with
# grep featurecount 3-RNA-Browser/*/config.yaml

# The run the 3-slurm.sh script

################################################################################
# Status check which genomes have completed what part

status <-
  todo.browser |>
  pull(dir) |>
  unique() |>
  map(function(i) {
    # i <- 'Acaryochloris.marina.MBIC11017_txid329726_SAMN02604308_unknown'
    tibble(
      dir = i,
      expression = file.path('3-RNA-Browser', i, 'analysis/41_counts.tsv') |>
        file.exists()
    )
  }) |>
  bind_rows()

################################################################################
################################################################################