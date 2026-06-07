# TE expression analysis

This document describes an example workflow for quantifying TE family expression from RNA-seq data using STAR and TEtranscripts.

## Purpose

The goal is to map RNA-seq reads to a soybean reference genome, generate a TE-compatible GTF file, quantify gene and TE expression, and convert raw counts to CPM values.

## Required tools

- FastQC
- fastp
- gffread
- STAR
- TEtranscripts / TEcount
- EDTA utility `split_overlap.pl`
- `makeTEgtf.pl`
- R or Python for count-table merging and CPM calculation

## Required input files

- `Wm82.fasta`: reference genome
- `Wm82.gff3`: gene annotation file
- `Wm82.mod.EDTA.bed`: TE BED file generated from EDTA annotation
- paired-end RNA-seq FASTQ files, for example `Soy01-A_1.fastq.gz` and `Soy01-A_2.fastq.gz`

## Step 1. Quality control

```bash
fastqc -t 8 \
  Soy01-A_1.fastq.gz \
  Soy01-A_2.fastq.gz
```

## Step 2. Read trimming

```bash
fastp \
  -i Soy01-A_1.fastq.gz \
  -I Soy01-A_2.fastq.gz \
  -o Soy01-A_clean_1.fq.gz \
  -O Soy01-A_clean_2.fq.gz \
  -5 -W 4 -M 20 -f 5 -e 20 -q 20 -l 50 -w 8
```

## Step 3. Convert GFF3 gene annotation to GTF

```bash
gffread Wm82.gff3 -T -o Wm82.gtf
```

## Step 4. Build a STAR genome index

```bash
STAR \
  --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir Wm82_STAR_index \
  --genomeFastaFiles Wm82.fasta \
  --sjdbGTFfile Wm82.gtf \
  --sjdbOverhang 149 \
  --genomeSAindexNbases 13
```

## Step 5. Map RNA-seq reads with STAR

Multimapped reads are retained because TE families often have multiple genomic copies.

```bash
STAR \
  --genomeDir Wm82_STAR_index \
  --readFilesIn Soy01-A_clean_1.fq.gz Soy01-A_clean_2.fq.gz \
  --readFilesCommand zcat \
  --runThreadN 8 \
  --genomeLoad NoSharedMemory \
  --outFilterMultimapNmax 1000 \
  --alignSJoverhangMin 8 \
  --alignSJDBoverhangMin 1 \
  --outFilterMismatchNmax 999 \
  --outFilterMismatchNoverReadLmax 0.04 \
  --alignIntronMin 20 \
  --alignIntronMax 1000000 \
  --alignMatesGapMax 1000000 \
  --outFilterScoreMinOverLread 0.5 \
  --outFilterMatchNminOverLread 0.5 \
  --outSAMheaderHD '@HD VN:1.4 SO:coordinate' \
  --outSAMunmapped Within \
  --outFilterType BySJout \
  --outSAMattributes NH HI AS NM MD \
  --outSAMtype BAM SortedByCoordinate \
  --sjdbScore 1 \
  --limitBAMsortRAM 30000000000 \
  --winAnchorMultimapNmax 1000 \
  --outFileNamePrefix Soy01-A_
```

The main BAM output will be:

```text
Soy01-A_Aligned.sortedByCoord.out.bam
```

## Step 6. Prepare a TE GTF file for TEtranscripts

Low-complexity, simple-repeat, and satellite records are removed. TE identifiers are also cleaned to avoid special characters in the GTF attributes.

```bash
# Remove non-TE repeat categories and clean TE names/orientation fields
awk -F'\t' 'BEGIN{OFS="\t"}
$5 !~ /Low_complexity|Simple_repeat|Satellite/ {
  gsub("[.?]", "C", $9);
  gsub(":", "_", $4);
  gsub("\\.\\.", "_", $4);
  print;
}' Wm82.mod.EDTA.bed \
  > Wm82.TE.cleaned.bed

# Split overlapping TE entries before GTF conversion
perl split_overlap.pl \
  Wm82.TE.cleaned.bed \
  Wm82.TE.cleaned.split.bed

# Convert TE BED to TE GTF
perl makeTEgtf.pl \
  -c 1 -s 2 -e 3 -o 9 -n RepeatMasker -t 4 \
  -1 Wm82.TE.cleaned.split.bed \
  > Wm82.TE.gtf
```

## Step 7. Quantify gene and TE expression

```bash
TEcount \
  -b Soy01-A_Aligned.sortedByCoord.out.bam \
  --GTF Wm82.gtf \
  --TE Wm82.TE.gtf \
  --mode multi \
  --sortByPos \
  --outdir TEcount_output \
  --project Soy01-A
```

## Step 8. Convert raw counts to CPM

The following R code converts a raw count table to CPM. The first column is assumed to contain gene or TE IDs, and the remaining columns are sample counts.

```r
library(tidyverse)

counts <- read_tsv("merged_gene_TE_counts.tsv", show_col_types = FALSE)

id_col <- counts[[1]]
count_mat <- counts[, -1] %>% as.data.frame()

cpm_mat <- sweep(count_mat, 2, colSums(count_mat), FUN = "/") * 1e6
cpm <- bind_cols(ID = id_col, as_tibble(cpm_mat))

write_tsv(cpm, "merged_gene_TE_CPM.tsv")
```

## Step 9. Define transcriptionally expressed TE families

A TE family can be considered transcriptionally expressed if its CPM value is greater than 1 in at least one sample.

```r
library(tidyverse)

cpm <- read_tsv("merged_gene_TE_CPM.tsv", show_col_types = FALSE)

expressed_te <- cpm %>%
  filter(if_any(-1, ~ .x > 1))

write_tsv(expressed_te, "expressed_TE_families_CPM_gt1.tsv")
```
