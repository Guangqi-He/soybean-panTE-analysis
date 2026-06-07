# TE distance to TSS and gene annotation

This document describes an example workflow for calculating the distance from each TE to the nearest transcription start site (TSS) and annotating TE positions relative to genes.

## Purpose

The goal is to determine whether TE copies are located upstream, downstream, within, or near genes, and to obtain signed distances between TE copies and the nearest TSS.

## Required tools

- BEDtools
- VEP
- bgzip and tabix
- samtools
- awk, grep, sort

## Required input files

- `Wm82.gff3`: gene annotation file
- `Wm82.fasta`: genome sequence
- `Wm82.TE.split.bed`: non-overlapping TE BED file

## Step 1. Generate a TSS BED file from gene annotation

The TSS is defined as the gene start coordinate for genes on the `+` strand and the gene end coordinate for genes on the `-` strand.

```bash
awk 'BEGIN{OFS="\t"}
$3 == "gene" {
  if ($7 == "+") {
    print $1, $4 - 1, $4, $9, ".", $7
  } else if ($7 == "-") {
    print $1, $5 - 1, $5, $9, ".", $7
  }
}' Wm82.gff3 \
  | sort -k1,1 -k2,2n \
  > Wm82.gene.TSS.bed
```

## Step 2. Prepare a sorted TE BED file

This example assumes that the TE BED file contains at least five columns: chromosome, start, end, TE ID, and TE superfamily.

```bash
grep -v -F '(' Wm82.TE.split.bed \
  | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5}' \
  | sort -k1,1 -k2,2n \
  > Wm82.TE.sorted.bed
```

## Step 3. Calculate signed distance to the nearest TSS

`bedtools closest -D b` reports signed distances relative to the feature in file `-b`. This makes it possible to distinguish upstream and downstream positions relative to the TSS.

```bash
bedtools closest \
  -a Wm82.TE.sorted.bed \
  -b Wm82.gene.TSS.bed \
  -D b \
  > Wm82.TE.nearest_TSS.signed.tsv
```

## Step 4. Prepare TE records for VEP annotation

VEP can annotate a six-column Ensembl-style input file with the following columns:

```text
chromosome  start  end  allele  strand  identifier
```

For TE interval annotation, the allele field is only used as a placeholder.

```bash
awk 'BEGIN{OFS="\t"}
{
  strand = "+";
  print $1, $2, $3, "DEL", strand, $4
}' Wm82.TE.sorted.bed \
  > Wm82.TE.vep.input.txt
```

## Step 5. Prepare GFF3 and FASTA files for VEP

The GFF3 file should be sorted, compressed, and indexed. The FASTA file should also be indexed.

```bash
grep -v '^#' Wm82.gff3 \
  | sort -k1,1 -k4,4n -k5,5n \
  | bgzip -c \
  > Wm82.gff3.gz

tabix -p gff Wm82.gff3.gz

bgzip -c Wm82.fasta > Wm82.fasta.gz
samtools faidx Wm82.fasta.gz
```

## Step 6. Run VEP annotation

```bash
vep \
  --format ensembl \
  -i Wm82.TE.vep.input.txt \
  --gff Wm82.gff3.gz \
  --fasta Wm82.fasta.gz \
  --stats_file Wm82.TE.vep.stats \
  --output_file Wm82.TE.vep.annotation.txt
```

## Step 7. Remove redundant transcript-level annotations

A TE may overlap multiple transcripts from the same gene. The following command keeps unique combinations of TE ID, genomic position, gene ID, and consequence.

```bash
grep -v '^#' Wm82.TE.vep.annotation.txt \
  | cut -f1,2,4,7 \
  | grep -v -F ')n' \
  | sort -u \
  > Wm82.TE.vep.annotation.unique.tsv
```

## Notes

- Use non-overlapping TE intervals before converting TE annotations to VEP input.
- The `--format ensembl` option is required for six-column VEP input.
- For signed TSS distances, the input BED files must be sorted by chromosome and start coordinate.

