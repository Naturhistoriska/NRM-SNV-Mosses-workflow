#!/usr/bin/env -S Rscript --vanilla

# Description: Script for PCA and basic clustering of populations
# By: Jason Hill, modifications by JN
# Last modified: ons apr 26, 2023  12:58

# Input files
sampleinfo <- "data/sample-population.tsv"
recodevcf <- "data/output_all.MAC1.recode.vcf"

# Output files:
treeplot <- "img/clustering.png"
pcaplot1 <- "img/pca_12.png"
pcaplot2 <- "img/pca_all.png"

# Settings
mycols <- c("red", "darkgreen", "blue")
gdsfile <- "output_all.MAC1.recode.gds"
set.seed(100)

# Packages
library(gdsfmt)
library(SNPRelate)
library(tidyverse)
library(ggtree)
library(ape)

# Read sample info
samples.info <- read_tsv(
    sampleinfo,
    col_names = c("sample", "population"))

# VCF to GDS
snpgdsVCF2GDS(
    recodevcf,
    gdsfile,
    method = "biallelic.only")

snps_MAC1.gds <- snpgdsOpen(gdsfile)

# Hierachical clustering
snps_MAC1.ibs <- snpgdsHCluster(
    snpgdsIBS(snps_MAC1.gds,num.thread = 2,
        autosome.only = FALSE,
        missing.rate = 0.5,
        maf = 0.1))

snps_MAC1.rv <- snpgdsCutTree(snps_MAC1.ibs)

tree <- ggtree(as.phylo(as.hclust(snps_MAC1.rv$dendrogram))) +
      geom_treescale() +
      xlim(0, 0.1)

png(treeplot)

tree %<+% samples.info +
    geom_tiplab(aes(fill = factor(population)),
        color = "black",
        geom = "label",
        size = 2)

dev.off()

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

png(pcaplot1)

plot(snps_MAC1.pop_pca$EV2,
     snps_MAC1.pop_pca$EV1,
     col = mycols[as.integer(snps_MAC1.pop_pca$pop)],
     xlab = "eigenvector 2",
     ylab = "eigenvector 1")

legend("topright",
       legend = levels(snps_MAC1.pop_pca$pop),
       pch = "o",
       col = mycols,
       title = "Population")

dev.off()

png(pcaplot2)

plot(snps_MAC1.pca,
    1:4,
    col = as.integer(snps_MAC1.pop_pca$pop))

dev.off()

snpgdsClose(snps_MAC1.gds)

unlink(gdsfile)

