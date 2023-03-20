# NRM SNV Mosses workflow

- Last modified: mån mar 20, 2023  12:42
- Sign: Johan Nylander

## Fastq filtering and Mapping

Set up directories

    $ mkdir -p /home/nylander/run/snv/eager/{data,reference}
    $ cd /home/nylander/run/snv/eager

Reference (genome)

    $ cp /path/to/reference/genome /home/nylander/run/snv/eager/reference/ref.fasta

Input data as fastq.gz files

    $ cd /home/nylander/run/snv/eager/data
    $ find /home/nylander/run/snv/data/P27213 -name '*_001.fastq.gz' -exec ln -s {} . \;
    $ cd /home/nylander/run/snv/eager

Eager parameters (different from default) in file `nf-params.json`:

    $ cat nf-params.json
    {
        "input": "data/*_{R1,R2}_*.fastq.gz",
        "fasta": "reference/ref.fasta",
        "fasta_index": "reference/ref.fasta.fai",
        "skip_damage_calculation": true,
        "mergedonly": true,
        "mapper": "bwamem",
        "run_genotyping": true,
        "genotyping_tool": "freebayes",
        "freebayes_p": 1
    }

Start eager pipeline

    $ screen -S eager
    $ cd /home/nylander/run/snv/eager
    $ nextflow run nf-core-eager-2.4.6/workflow \
      -name run_2 \
      -profile singularity \
      -params-file nf-params.json

Clean up

    $ rm -r work
    $ nextflow clean -f -k

---

## Haplotype calling

Set up directories

    $ mkdir -p /home/nylander/run/snv/freebayes/reference

Reference

    $ ln -s /home/nylander/run/snv/eager/reference/ref.fasta /home/nylander/run/snv/freebayes/reference/ref.fasta

List of bam files (should have indexes)

    $ cd /home/nylander/run/snv/freebayes
    $ find /home/nylander/run/snv/eager/results -name '*_rmdup.bam' > bam.list

Run freebayes, followed by sorting the VCF

    $ cd /home/nylander/run/snv/freebayes
    $ freebayes \
        --fasta-reference reference/ref.fasta \
        --targets cov30.bed \
        --bam-list bamlist.txt \
        --ploidy 1 \
        --gvcf \
        --skip-coverage 700 \
        --vcf freebayes_gvcf.vcfs

    $ bcftools sort freebayes_gvcf.vcf -o freebayes_gvcf_sorted.vcf.gz -Oz

---

## Filter VCF (Jason Hill)

The VCF file derived from the Freebayes join variant calling run contained 4039951 variant sites. Inspection showed that a near majority of these sites were bizare looking indels, however most of these had an allele frequency of 0 so they must have been included by freebayes for some reason even though none of the samples carried the variant indel. Some basic filtering of the VCF produced a more reasonable set.

    $ vcftools --gzvcf freebayes_gvcf_sorted.vcf.gz --minDP 20 --minQ 20 --recode-INFO-all --out filtered_GQ20_DP20

Require genotype quality > 20 and depth > 20. The depth filter was redundant since only sites with depth > 30 were supplied anyway.

    $ vcftools --vcf filtered_GQ20_DP20.vcf --remove-indels --recode --recode-INFO-all --out output_snps-only.vcf

How many SNPs? After filtering, kept 797981 out of a possible 1501403 Sites

    $ vcftools --vcf filtered_GQ20_DP20.vcf --keep-only-indels --recode --recode-INFO-all --out output_indels-only.vcf

How many indels? After filtering, kept 703422 out of a possible 1501403 Sites

    $ vcftools --vcf filtered_GQ20_DP20.vcf --recode --recode-INFO-all --mac 1 --out output_all.MAC1.vcf

We'll use both SNPs and indels for now, but require that at least one sample actually has the alternate allele.After filtering, kept 398100 out of a possible 1501403 Sites. There was a very large reduction in the number of variant sites, but those that remain look much more reasonable. This averages out to one variant every 550bp, which seems somewhat low? We’ll use this file going forward

## Clustering of samples and search for structure (Jason Hill)

Infile `sample-population.tsv`:

    $ cat sample-population.tsv
    P27213_101_S53_L002	popA
    P27213_102_S54_L002	popA
    P27213_103_S55_L002	popA
    P27213_104_S56_L002	popA
    P27213_105_S57_L002	popA
    P27213_106_S58_L002	popA
    P27213_107_S59_L002	popA
    P27213_108_S60_L002	popA
    P27213_109_S61_L002	popA
    P27213_110_S62_L002	popA
    P27213_111_S63_L002	popB
    P27213_112_S64_L002	popB
    P27213_113_S65_L002	popB
    P27213_114_S66_L002	popB
    P27213_115_S67_L002	popB
    P27213_116_S68_L002	popB
    P27213_117_S69_L002	popB
    P27213_118_S70_L002	popB
    P27213_119_S71_L002	popB
    P27213_120_S72_L002	popB
    P27213_121_S73_L002	popC
    P27213_122_S74_L002	popC
    P27213_123_S75_L002	popC
    P27213_124_S76_L002	popC
    P27213_125_S77_L002	popC
    P27213_126_S78_L002	popC
    P27213_127_S79_L002	popC
    P27213_128_S80_L002	popC
    P27213_129_S81_L002	popC
    P27213_130_S82_L002	popC

Commands in R

```R

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

```

## Figures

---

![Hierarchical clustering](img/clustering.png)

---

![PCA 1, 2](img/pca12.png)

---

![PCA all](img/pcaall.png)

---

## Ackowledgements

Thanks to Jason Hill [NBIS.se](https://nbis.se/about/staff/jason-hill/) for advice on analyses.

