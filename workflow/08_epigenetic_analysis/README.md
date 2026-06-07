# Epigenetic analysis

This document describes an example workflow for integrating TE-derived variants with chromatin accessibility, histone modification, DNA methylation, and gene expression data.

## Purpose

The goal is to construct a mini linear pangenome, map epigenomic and transcriptomic reads to it, and evaluate whether TE-derived variant regions are associated with chromatin accessibility, DNA methylation, histone marks, or gene expression changes.

## Required tools

- meryl
- winnowmap
- samtools
- minipileup
- bcftools
- liftoff
- EDTA
- bowtie2
- sambamba
- BEDtools
- MACS3
- IDR
- Bismark
- STAR
- gffread
- R

## Required input files

- `Wm82.fasta`: reference genome
- `glysoW05.fasta`: query genome used to construct a mini linear pangenome
- `Wm82.gff3`: gene annotation of the reference genome
- `Wm82.cds.rep.fa`: representative CDS sequence file
- `panTE_DeepTE_final.lib.fa`: final pan-genome TE library
- paired-end ATAC-seq, ChIP-seq, BS-seq, and RNA-seq FASTQ files

## Step 1. Call W05 insertions relative to Wm82

```bash
ref=Wm82.fasta
query=glysoW05.fasta
sample=glysoW05

meryl count k=19 output Wm82.meryl ${ref}
meryl print greater-than distinct=0.9998 Wm82.meryl > Wm82_k19.txt

winnowmap -ax asm5 --eqx -t 24 -W Wm82_k19.txt \
  ${ref} ${query} \
  > ${sample}.sam

samtools sort -O BAM -o ${sample}.bam ${sample}.sam
samtools index ${sample}.bam

minipileup -f ${ref} -v -c -C -q 0 -s 0 -a 0 -l 1000 \
  ${sample}.bam \
  > ${sample}.PAV.vcf

bcftools norm -m -any -f ${ref} ${sample}.PAV.vcf \
  -Oz -o ${sample}.PAV.norm.vcf.gz

bcftools index -t ${sample}.PAV.norm.vcf.gz
```

## Step 2. Keep insertion alleles longer than 100 bp

```bash
bcftools view -H ${sample}.PAV.norm.vcf.gz \
  | awk 'BEGIN{FS=OFS="\t"}
{
  ref = $4;
  alt = $5;
  n = split(alt, a, ",");
  keep = 0;
  for (i = 1; i <= n; i++) {
    if (length(a[i]) > length(ref) && length(a[i]) - length(ref) >= 100) {
      keep = 1;
      break;
    }
  }
  if (keep == 1) print;
}' > ${sample}.INS100.body.vcf

bcftools view -h ${sample}.PAV.norm.vcf.gz > ${sample}.INS100.header.vcf
cat ${sample}.INS100.header.vcf ${sample}.INS100.body.vcf \
  | grep -v 'NNNNNNNNNN' \
  > ${sample}.INS100.noN.vcf

bgzip -f ${sample}.INS100.noN.vcf
bcftools index -t ${sample}.INS100.noN.vcf.gz
```

## Step 3. Construct a mini linear pangenome

```bash
bcftools consensus \
  -f ${ref} \
  -H 1 \
  -o minilinerpan.fa \
  ${sample}.INS100.noN.vcf.gz

seqkit stat minilinerpan.fa
```

## Step 4. Lift reference gene annotation to the mini linear pangenome

```bash
liftoff \
  -g Wm82.gff3 \
  -o minilinerpan.gff3 \
  -p 24 \
  -polish \
  -copies \
  minilinerpan.fa \
  ${ref}
```

## Step 5. Annotate TEs in the mini linear pangenome

```bash
EDTA.pl \
  --genome minilinerpan.fa \
  --overwrite 0 \
  --sensitive 1 \
  --anno 1 \
  --evaluate 0 \
  --force 1 \
  --curatedlib panTE_DeepTE_final.lib.fa \
  --cds Wm82.cds.rep.fa \
  --threads 32
```

## Step 6. Identify W05- and W82-specific insertion regions on the mini linear pangenome

After constructing the mini linear pangenome, align both Wm82 and W05 back to the mini linear pangenome. A deletion in Wm82 relative to the mini pangenome corresponds to a W05 insertion, whereas a deletion in W05 corresponds to a Wm82 insertion.

```bash
meryl count k=19 output minipan.meryl minilinerpan.fa
meryl print greater-than distinct=0.9998 minipan.meryl > minipan_k19.txt

winnowmap -ax asm5 --eqx -t 24 -W minipan_k19.txt \
  minilinerpan.fa Wm82.fasta \
  > Wm82.minipan.sam

winnowmap -ax asm5 --eqx -t 24 -W minipan_k19.txt \
  minilinerpan.fa glysoW05.fasta \
  > glysoW05.minipan.sam

samtools sort -O BAM -o Wm82.minipan.bam Wm82.minipan.sam
samtools sort -O BAM -o glysoW05.minipan.bam glysoW05.minipan.sam

minipileup -f minilinerpan.fa -v -c -C -q 0 -s 0 -a 0 -l 1000 \
  Wm82.minipan.bam \
  > Wm82.minipan.PAV.vcf

minipileup -f minilinerpan.fa -v -c -C -q 0 -s 0 -a 0 -l 1000 \
  glysoW05.minipan.bam \
  > glysoW05.minipan.PAV.vcf
```

Extract deletion intervals from each VCF. Only deletions at least 50 bp are retained.

```bash
awk 'BEGIN{FS=OFS="\t"}
/^#/ {next}
length($4) > length($5) {
  len = length($4) - length($5);
  start = $2;
  end = $2 + len;
  if (len >= 50) print $1, start, end, len;
}' Wm82.minipan.PAV.vcf > W05_insertion_regions.bed

awk 'BEGIN{FS=OFS="\t"}
/^#/ {next}
length($4) > length($5) {
  len = length($4) - length($5);
  start = $2;
  end = $2 + len;
  if (len >= 50) print $1, start, end, len;
}' glysoW05.minipan.PAV.vcf > W82_insertion_regions.bed
```

## Step 7. Map ATAC-seq or ChIP-seq reads to the mini linear pangenome

The ATAC-seq example uses local alignment. ChIP-seq can use `--very-sensitive` instead of `--very-sensitive-local`.

```bash
bowtie2-build --threads 24 minilinerpan.fa minipan

bowtie2 -x minipan \
  -1 W05_ATAC_rep1_1.fq.gz \
  -2 W05_ATAC_rep1_2.fq.gz \
  -X 2000 --dovetail --no-mixed --no-discordant \
  --very-sensitive-local -p 32 --no-unal \
  | samtools view -b -q 30 -F 1804 -@ 4 \
  | samtools sort -O BAM -@ 4 -o W05_ATAC_rep1.sort.bam

sambamba markdup --remove-duplicates --nthreads 8 \
  W05_ATAC_rep1.sort.bam \
  W05_ATAC_rep1.rmDup.bam
```

## Step 8. Call ATAC-seq peaks

```bash
bedtools bamtobed -i W05_ATAC_rep1.rmDup.bam > W05_ATAC_rep1.bed

macs3 callpeak \
  -g 1028537865 \
  -t W05_ATAC_rep1.bed \
  --name W05_ATAC_rep1 \
  --format BED \
  --keep-dup all \
  --qvalue 0.05 \
  --cutoff-analysis \
  --broad \
  --nomodel --shift -100 --extsize 200
```

To obtain reproducible broad peaks from two replicates:

```bash
awk '$9 >= 1.3' W05_ATAC_rep1_peaks.broadPeak > W05_ATAC_rep1.q05.broadPeak
awk '$9 >= 1.3' W05_ATAC_rep2_peaks.broadPeak > W05_ATAC_rep2.q05.broadPeak

bedtools intersect -f 0.5 -r \
  -a W05_ATAC_rep1.q05.broadPeak \
  -b W05_ATAC_rep2.q05.broadPeak \
  -wa \
  > W05_ATAC_reproducible_peaks.bed
```

## Step 9. Call ChIP-seq peaks

Example for a narrow histone mark:

```bash
macs3 callpeak \
  -g 1028537865 \
  -t W05_H3K27ac_rep1.rmDup.bam \
  -c W05_input_rep1.rmDup.bam W05_input_rep2.rmDup.bam \
  --name W05_H3K27ac_rep1 \
  --format BAMPE \
  --keep-dup all \
  --qvalue 0.05 \
  --cutoff-analysis
```

Example for a broad histone mark such as H3K27me3:

```bash
macs3 callpeak \
  --broad \
  -g 1028537865 \
  -t W05_H3K27me3_rep1.rmDup.bam \
  -c W05_input_rep1.rmDup.bam W05_input_rep2.rmDup.bam \
  --name W05_H3K27me3_rep1 \
  --format BAMPE \
  --keep-dup all \
  --qvalue 0.05 \
  --cutoff-analysis
```

## Step 10. Map BS-seq reads and extract DNA methylation information

```bash
mkdir methyindex
ln -s minilinerpan.fa methyindex/

bismark_genome_preparation \
  --bowtie2 \
  --parallel 16 \
  methyindex

bismark \
  -N 0 -L 20 --un --ambiguous --bam --bowtie2 -p 8 \
  -o BS_output \
  --fastq --genome methyindex \
  -1 W05_BS_rep1_1.fq.gz \
  -2 W05_BS_rep1_2.fq.gz

deduplicate_bismark \
  -bam -p \
  --output_dir BS_output \
  BS_output/W05_BS_rep1_1_bismark_bt2_pe.bam

bismark_methylation_extractor \
  -p --parallel 4 --gzip --bedGraph --no_overlap \
  --comprehensive --counts --CX_context --cytosine_report \
  --buffer_size 20G \
  --genome_folder methyindex \
  -o BS_output/W05_BS_rep1_report \
  BS_output/W05_BS_rep1_1_bismark_bt2_pe.deduplicated.bam
```

## Step 11. Summarize methylation counts over TE-derived variant regions

The Bismark `.cov.gz` file is converted to a five-column BED-like file: chromosome, start, end, methylated read count, and unmethylated read count.

```bash
zcat W05_BS_rep1.bismark.cov.gz \
  | awk 'BEGIN{OFS="\t"} {print $1, $2 - 1, $3, $5, $6}' \
  > W05_BS_rep1.mCuC.bed

bedtools map \
  -a W05_insertion_regions.bed \
  -b W05_BS_rep1.mCuC.bed \
  -c 4,5 \
  -o sum,sum \
  > W05_BS_rep1.W05_insertion_region_counts.bed
```

## Step 12. Map RNA-seq reads and quantify gene expression

```bash
gffread minilinerpan.gff3 -T -o minilinerpan.gtf

STAR \
  --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir minipan_STAR_index \
  --genomeFastaFiles minilinerpan.fa \
  --sjdbGTFfile minilinerpan.gtf \
  --sjdbOverhang 149 \
  --genomeSAindexNbases 13

STAR \
  --genomeDir minipan_STAR_index \
  --runThreadN 8 \
  --readFilesIn W05_RNA_rep1_1.fq.gz W05_RNA_rep1_2.fq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix W05_RNA_rep1 \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMstrandField intronMotif \
  --outSAMattributes All \
  --quantMode GeneCounts \
  --outFilterScoreMinOverLread 0.5 \
  --outFilterMatchNminOverLread 0.5 \
  --outFilterIntronMotifs RemoveNoncanonical
```
