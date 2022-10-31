library(tidyverse)

pacbio_file <- "./result/annotations/pacbio/jtk/pacbio.tsv"
pacbio <- read_tsv(pacbio_file, col_names = FALSE)
ler_file <- "./result/annotations/ler/jtk/ler.tsv"
ler <- read_tsv(ler_file, col_names = FALSE)

merged_data <- full_join(
    pacbio %>% select(X2, X3) %>% mutate(PacBio = TRUE),
    ler %>% select(X2, X3) %>% mutate(Ler = TRUE),
    by = c("X2", "X3"),
)
outfile <- "./result/analysis/pacbio_vs_ler.tsv"
write_tsv(merged_data, outfile)
