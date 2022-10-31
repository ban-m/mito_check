#!/bin/env R
library(tidyverse)
loadNamespace("cowplot")
kmer_counts <- read_tsv("./result/repeat_annotations/kmer_counts.tsv")
g <- kmer_counts %>% filter(OccInGenome < 50) %>% ggplot() +
    geom_point(aes(x = OccInGenome, y = NumOfKmer)) +
    facet_wrap(vars(K)) +
    cowplot::theme_cowplot() +
    labs(x = "Occurence in the genomes", y = "Number of Kmer")
cowplot::ggsave2(filename = "./pics/kmer_histogram.pdf", plot = g)
cowplot::ggsave2(filename = "./pics/kmer_histogram.png", plot = g)

file <- "./result/repeat_annotations/dist_to_nearest_repeats.tsv"
dist_to_repeats <- read_tsv(file)
g <- dist_to_repeats %>% ggplot() +
    geom_histogram(aes(x = Distance, y = ..density..), bins = 60) +
    facet_grid(IsNull ~ .) +
    cowplot::theme_cowplot() +
    labs(x = "Distance from the nearest repeats (bp)", y = "Frequency")
cowplot::ggsave2(filename = "./pics/dist_to_repeats.pdf", plot = g)
cowplot::ggsave2(filename = "./pics/dist_to_repeats.png", plot = g)

file <- "./result/repeat_annotations/average_dist_to_repeats.tsv"
average_dist <- read_tsv(file)
obs_ave <- average_dist %>% filter(IsNull == "Obs") %>% pull(Distance)
g <- average_dist %>% filter(IsNull == "Sampled") %>% ggplot() +
    geom_histogram(aes(x = Distance), bins = 60) +
    geom_vline(xintercept = obs_ave) +
    cowplot::theme_cowplot() +
    labs(x = "Average distance from the nearest repeats(bp)", y = "Frequency")
cowplot::ggsave2(filename = "./pics/average_to_repeats.pdf", plot = g)
cowplot::ggsave2(filename = "./pics/average_to_repeats.png", plot = g)