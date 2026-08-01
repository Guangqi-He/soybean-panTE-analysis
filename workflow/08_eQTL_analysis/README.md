# TAP/TIP eQTL Analysis

This workflow genotypes TE presence/absence variants from resequencing data, processes RNA-seq expression data, and performs cis- and trans-eQTL mapping with TensorQTL.

## Required input

- Paired-end resequencing FASTQ files
- Paired-end RNA-seq FASTQ files
- PanGenie index and the corresponding variant callset
- STAR genome index
- Gene coordinates in BED format
- Covariate matrix
- Sample name mapping file

Main software: fastp, PanGenie, BCFtools, PLINK2, STAR, DESeq2, TensorQTL, pandas, and pyarrow.

Sample lists contain one sample ID per line. FASTQ files are assumed to follow this naming scheme:

```text
sample_f1.fq.gz
sample_r2.fq.gz
```

## 1. TAP/TIP genotyping

### 1.1 Read trimming

```bash
while read -r id; do fastp --in1 ${id}_f1.fq.gz --in2 ${id}_r2.fq.gz --out1 01.qc/${id}_1.fq.gz --out2 01.qc/${id}_2.fq.gz --trim_front1 5 --trim_front2 5; done < resequencing.list
```

### 1.2 PanGenie genotyping

```bash
while read -r id; do PanGenie -f pangenie_index/panel -i <(zcat 01.qc/${id}_1.fq.gz 01.qc/${id}_2.fq.gz) -s ${id} -o pangenie_out/${id} -j 24 -t 24; bgzip -f pangenie_out/${id}_genotyping.vcf; tabix -f -p vcf pangenie_out/${id}_genotyping.vcf.gz; done < resequencing.list
```

Convert each PanGenie result to the biallelic callset representation and mask low-quality genotypes:

```bash
while read -r id; do zcat pangenie_out/${id}_genotyping.vcf.gz | python3 convert-to-biallelic.py panel_callset.vcf.gz | bgzip -c > pangenie_out/${id}.callset.vcf.gz; bcftools +setGT pangenie_out/${id}.callset.vcf.gz -Oz -o pangenie_out/${id}.callset.gq200.vcf.gz -- -t q -n . -i 'FMT/GQ<200'; tabix -f -p vcf pangenie_out/${id}.callset.gq200.vcf.gz; done < resequencing.list
```

Merge all samples:

```bash
awk '{print "pangenie_out/"$1".callset.gq200.vcf.gz"}' resequencing.list > vcf.list
bcftools merge -l vcf.list --threads 4 -Oz -o pangenie_all.vcf.gz
bcftools index -t pangenie_all.vcf.gz
```

The `GQ >= 200` cutoff was used in this analysis and can be adjusted for other datasets.

## 2. RNA-seq processing

### 2.1 Read trimming

```bash
while read -r id; do fastp --in1 ${id}_f1.fq.gz --in2 ${id}_r2.fq.gz --out1 01.qc/${id}_1.fq.gz --out2 01.qc/${id}_2.fq.gz --trim_front1 5 --trim_front2 5; done < RNA.list
```

### 2.2 Alignment and gene counting

```bash
while read -r id; do STAR --runThreadN 6 --genomeDir STAR_index --readFilesIn 01.qc/${id}_1.fq.gz 01.qc/${id}_2.fq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix 02.map/${id}. --quantMode GeneCounts; done < RNA.list
```

Merge the `ReadsPerGene.out.tab` files into one count matrix. STAR column 2 is used for unstranded RNA-seq; use column 3 or 4 for stranded libraries when appropriate.

```bash
python3 merge_counts.py RNA.list 02.map 2 > gene_counts.tsv
python3 statSTARmappingResult.py RNA.list mappingsStat.tsv
```

### 2.3 Expression normalization

The raw count matrix is normalized with DESeq2, transformed with VST, and then converted by rank-based inverse normal transformation for eQTL mapping.

Save the following code as `normalize_expression.R`:

```r
library(DESeq2)

count_matrix <- read.table("gene_counts.tsv", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
count_matrix <- as.matrix(count_matrix)
storage.mode(count_matrix) <- "integer"

sample_map <- read.table("sample_map.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
new_names <- sample_map[[2]][match(colnames(count_matrix), sample_map[[1]])]
hit <- !is.na(new_names)
colnames(count_matrix)[hit] <- new_names[hit]
stopifnot(!anyDuplicated(colnames(count_matrix)))

coldata <- data.frame(row.names = colnames(count_matrix), group = factor(rep("all", ncol(count_matrix))))
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = coldata, design = ~1)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- estimateSizeFactors(dds)

normalized_counts <- counts(dds, normalized = TRUE)
vst_matrix <- assay(vst(dds, blind = TRUE))

inverse_normal_transform <- function(x) {
  r <- rank(x, ties.method = "average", na.last = "keep")
  qnorm((r - 3 / 8) / (sum(!is.na(x)) + 1 / 4))
}

expression_int <- t(apply(vst_matrix, 1, inverse_normal_transform))
rownames(expression_int) <- rownames(vst_matrix)
colnames(expression_int) <- colnames(vst_matrix)

write.table(normalized_counts, "expression.normalized.tsv", sep = "\t", quote = FALSE, col.names = NA)
write.table(vst_matrix, "expression.vst.tsv", sep = "\t", quote = FALSE, col.names = NA)
write.table(expression_int, "expression.INT.tsv", sep = "\t", quote = FALSE, col.names = NA)
```

Run the script:

```bash
Rscript normalize_expression.R
```

## 3. Prepare TensorQTL input

### 3.1 Prepare the genotype file

Rename samples, sort the VCF, and assign a unique ID to each variant:

```bash
bcftools reheader -s sample_map.txt -o genotypes.renamed.vcf.gz pangenie_all.vcf.gz
bcftools sort genotypes.renamed.vcf.gz -Oz -o genotypes.sorted.vcf.gz
bcftools view -Ov genotypes.sorted.vcf.gz | awk 'BEGIN{FS=OFS="\t"} /^#/{print; next} {$3=$1":"$2":"++n; print}' | bgzip -c > genotypes.vcf.gz
tabix -f -p vcf genotypes.vcf.gz
bcftools query -l genotypes.vcf.gz > genotype_samples.list
```

### 3.2 Prepare phenotype and covariate files

The phenotype BED and covariate matrix must contain the same samples as the genotype file and in the same order.

```bash
python3 make_tensorqtl_inputs.py --gene-bed genes.bed --expr expression.INT.tsv --cov covariates.txt --vcf-samples genotype_samples.list --out-bed expression.tensorqtl.bed.gz --out-cov covariates.tensorqtl.txt
```

Convert the VCF to PLINK2 format:

```bash
plink2 --vcf genotypes.vcf.gz --make-pgen --out genotypes --allow-extra-chr
```

## 4. eQTL mapping

Create output directories:

```bash
mkdir -p cis cis_independent cis_nominal trans
```

Run cis-eQTL mapping with a 1 Mb window, 10,000 permutations, and a minor allele frequency cutoff of 0.01:

```bash
python3 -m tensorqtl genotypes expression.tensorqtl.bed.gz sv_eqtl --covariates covariates.tensorqtl.txt --mode cis --window 1000000 --permutations 10000 --chunk_size chr --maf_threshold 0.01 -o cis
```

Identify conditionally independent cis-eQTLs:

```bash
python3 -m tensorqtl genotypes expression.tensorqtl.bed.gz sv_eqtl --covariates covariates.tensorqtl.txt --cis_output cis/sv_eqtl.cis_qtl.txt.gz --mode cis_independent -o cis_independent
```

Run nominal cis-eQTL mapping:

```bash
python3 -m tensorqtl genotypes expression.tensorqtl.bed.gz sv_eqtl --covariates covariates.tensorqtl.txt --mode cis_nominal --window 1000000 --chunk_size chr -o cis_nominal
```

Run trans-eQTL mapping:

```bash
python3 -m tensorqtl genotypes expression.tensorqtl.bed.gz sv_eqtl --covariates covariates.tensorqtl.txt --mode trans --pval_threshold 1e-5 --batch_size 20000 --chunk_size 5000 --output_text -o trans
```

## 5. Output processing

TensorQTL nominal cis results are stored as Parquet files. Convert them to compressed text files when needed:

```python
import glob
import os
import pandas as pd

for infile in glob.glob("cis_nominal/*.parquet"):
    outfile = os.path.splitext(infile)[0] + ".txt.gz"
    pd.read_parquet(infile).to_csv(outfile, sep="\t", index=False, compression="gzip")
```

A project-specific post-processing script can be used to combine significant cis, nominal cis, and trans associations:

```bash
python3 tensorqtl_postprocess.py --cis cis/sv_eqtl.cis_qtl.txt.gz --cis-nominal-dir cis_nominal --trans trans/sv_eqtl.trans_qtl_pairs.txt.gz -o tensorqtl_post --fdr 0.05 --chunksize 500000
```

## Notes

- Sample names must be identical across genotype, expression, and covariate files.
- The sample order in the phenotype and covariate files must match `genotype_samples.list`.
- Chromosome names must be consistent among the VCF, gene BED, and expression BED files.
- The gene IDs in `genes.bed` must match those in the expression matrix.
- Filtering thresholds should be adjusted according to sequencing depth, sample size, and study design.
