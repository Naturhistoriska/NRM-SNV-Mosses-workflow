# Description: Script for PCA and basic clustering of populations
# By: Jason Hill
# Version: tis 25 apr 2023 15:47:34
#
# Input:
#   - sample-population.tsv
#   - output_all.MAC1.vcf.recode.vcf

library(gdsfmt)
library(SNPRelate)
library(tidyverse)
library(ggtree)
library(ape)

samples.info <- read_tsv(
    "sample-population.tsv",
    col_names = c("sample", "population"))

# VCF to GDS
snpgdsVCF2GDS(
    "output_all.MAC1.vcf.recode.vcf",
    "output_all.MAC1.gds",
    method = "biallelic.only")

snps_MAC1.gds <- snpgdsOpen("output_all.MAC1.gds")

set.seed(100)

# Hierachical clustering
snps_MAC1.ibs <- snpgdsHCluster(
    snpgdsIBS(snps_MAC1.gds,num.thread = 2,
        autosome.only = FALSE,
        missing.rate = 0.5,
        maf = 0.1))

snps_MAC1.rv <- snpgdsCutTree(snps_MAC1.ibs)

p1 <- ggtree(as.phylo(
    as.hclust(snps_MAC1.rv$dendrogram))) +
    geom_treescale() +
    xlim(0, 0.1)

p1 %<+% samples.info +
    geom_tiplab(aes(fill = factor(population)),
        color = "black",
        geom = "label")

# PCA
snps_MAC1.pca <- snpgdsPCA(
    snps_MAC1.gds,
    autosome.only = FALSE,
    missing.rate = 0.5,
    maf = 0.1)

snps_MAC1.pop_pca <- data.frame(
    sample.id = snps_MAC1.pca$sample.id,
    pop = factor(samples.info$population)[match(snps_MAC1.pca$sample.id,
        samples.info$sample)],
            EV1 = snps_MAC1.pca$eigenvect[,1],
            EV2 = snps_MAC1.pca$eigenvect[,2],
            stringsAsFactors = FALSE)

plot(snps_MAC1.pop_pca$EV2,
     snps_MAC1.pop_pca$EV1,
     col = as.integer(snps_MAC1.pop_pca$pop),
     xlab = "eigenvector 2",
     ylab = "eigenvector 1")

legend("right",
       legend = levels(snps_MAC1.pop_pca$pop),
       pch = "o",
       col = 1:4)

plot(snps_MAC1.pca,
    1:4,
    col = as.integer(snps_MAC1.pop_pca$pop))

snpgdsClose(snps_MAC1.gds)

