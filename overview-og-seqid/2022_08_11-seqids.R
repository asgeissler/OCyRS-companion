library(tidyverse)

# preliimnary seqids foo
# vs
# A new genomic taxonomy system for the Synechococcus collective
# Vinícius W. Salazar,Diogo A. Tschoeke,Jean Swings,Carlos A. Cosenza,Marta Mattoso,Cristiane C. Thompson,Fabiano L. Thompson

og.sid <- read_tsv('2022_08_11-seqid.tsv')
ref <- read_tsv('2022_08_11-emi15173-sup-0004-tables2.tsv')


bind_rows(
  transmute(og.sid, x = 'This study', avg),
  transmute(ref, x = 'Salazar et al.', avg = `Mean AAI`)
) %>%
  ggplot(aes(avg, col = x)) +
  stat_ecdf(size = 2) +
  scale_color_manual(values = c( '#56B4E9', '#E69F00'), name = NULL) +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  xlab('Average pairwise sequence id') +
  ylab('Emp. cum. density') +
  theme_bw(18) +
  theme(legend.position = 'bottom')


ks.test(
  og.sid$avg, ref$`Mean AAI`,
  alternative = 'less'
)
#  p-value < 2.2e-16

ggsave('2022_08_11-fig.jpeg', width = 8, height = 6)
