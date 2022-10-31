library(tidyverse)
loadNamespace("cowplot")

strains <- c("an1", "c24", "col0", "cvi", "eri", "kyo", "ler", "pacbio", "sha")

break_scale <- 50

plot_coverage <- function(strain) {
    data_file <- paste0("./result/", strain, "/jtk/", strain, ".coverage.tsv")
    outfile <- paste0("./pics/coverages/", strain, ".png")
    dataset <- read_tsv(data_file)
    length <- dim(dataset)[1]
    dataset <- dataset %>% slice(5:(length - 5))
    max_coverage <- max(dataset %>% pull(coverage))
    upper <- ceiling((max_coverage * 1.2) / break_scale) * break_scale
    num_breaks <- upper / break_scale + 1
    print(upper)
    print(num_breaks)
    breaks <- seq(0, upper, length = num_breaks)
    g <- dataset %>% ggplot() +
        geom_line(aes(x = position, y = coverage)) + facet_grid(contig ~ .) +
        scale_y_continuous(breaks = breaks, limits = c(0, upper)) +
        cowplot::theme_cowplot()
    cowplot::ggsave2(filename = outfile, plot = g)
}

strains %>% sapply(FUN = plot_coverage)
