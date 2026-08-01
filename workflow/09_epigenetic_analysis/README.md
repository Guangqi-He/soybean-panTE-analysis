# Wm82–W05 Mini-Pangenome Multi-Omics Analysis

This workflow uses Wm82 as the backbone and adds long W05 insertions relative to Wm82 to construct a mini linear pangenome. ATAC-seq, ChIP-seq, BS-seq, and RNA-seq reads are then mapped to the same genome to examine the relationship between TE presence/absence variants, epigenetic changes, and gene expression.

## 1. Input files

```text
Wm82.fasta
Wm82.gff3
Wm82.cds.rep.fa
glysoW05.fasta
panTE_DeepTE_final.lib.fa
```

FASTQ files should be named as follows:

```text
sample_1.fq.gz
sample_2.fq.gz
```

Prepare `ATAC.list`, `ChIP.list`, `BS.list`, and `RNAseq.list`, with one sample name per line and without the FASTQ suffix.

Main software: meryl, Winnowmap, SAMtools, minipileup, BCFtools, Liftoff, EDTA, Bowtie2, MACS3, BEDTools, featureCounts, Bismark, STAR, and DESeq2.

---

## 2. Construct the mini linear pangenome

### 2.1 Identify W05 insertions relative to Wm82

```bash
ref=Wm82.fasta
query=glysoW05.fasta

meryl count k=19 output Wm82.meryl ${ref}
meryl print greater-than distinct=0.9998 Wm82.meryl > Wm82_k19.txt
winnowmap -ax asm5 --eqx -t 24 -W Wm82_k19.txt ${ref} ${query} | samtools sort -@ 8 -O BAM -o glysoW05.to_Wm82.bam
samtools index glysoW05.to_Wm82.bam
minipileup -f ${ref} -v -c -C -q 0 -s 0 -a 0 -l 1000 glysoW05.to_Wm82.bam > glysoW05.PAV.vcf
bcftools norm -m -any -f ${ref} -Oz -o glysoW05.PAV.norm.vcf.gz glysoW05.PAV.vcf
bcftools index -t glysoW05.PAV.norm.vcf.gz
```

Keep insertions at least 100 bp long and remove alleles containing `N`:

```bash
( bcftools view -h glysoW05.PAV.norm.vcf.gz; bcftools view -H glysoW05.PAV.norm.vcf.gz | awk 'BEGIN{FS=OFS="\t"} $5!="." && $5!="*" && $5!~/^</ && $4!~/[Nn]/ && $5!~/[Nn]/ && length($5)-length($4)>=100' ) | bgzip -c > glysoW05.INS100.vcf.gz
bcftools index -t glysoW05.INS100.vcf.gz
```

### 2.2 Build the mini-pangenome

```bash
bcftools consensus -f ${ref} -s - -c Wm82_to_minipan.chain -o minilinerpan.fa glysoW05.INS100.vcf.gz
samtools faidx minilinerpan.fa
liftoff -g Wm82.gff3 -o minilinerpan.gff3 -p 24 -polish -copies minilinerpan.fa ${ref}
EDTA.pl --genome minilinerpan.fa --overwrite 0 --sensitive 1 --anno 1 --evaluate 0 --curatedlib panTE_DeepTE_final.lib.fa --cds Wm82.cds.rep.fa --threads 32
```

Convert the EDTA annotation to BED. Adjust the input filename if necessary:

```bash
awk 'BEGIN{FS=OFS="\t"} $0!~/^#/ {print $1,$4-1,$5,$3,$9}' minilinerpan.fa.mod.EDTA.TEanno.gff3 | sort -k1,1 -k2,2n > minilinerpan.TE.bed
```

---

## 3. Identify genome-specific sequences

Map Wm82 and W05 back to the mini-pangenome. Deletions in the Wm82 alignment represent W05-specific insertions, whereas deletions in the W05 alignment represent W82-specific sequences.

```bash
meryl count k=19 output minipan.meryl minilinerpan.fa
meryl print greater-than distinct=0.9998 minipan.meryl > minipan_k19.txt
winnowmap -ax asm5 --eqx -t 24 -W minipan_k19.txt minilinerpan.fa Wm82.fasta | samtools sort -@ 8 -O BAM -o Wm82.to_minipan.bam
winnowmap -ax asm5 --eqx -t 24 -W minipan_k19.txt minilinerpan.fa glysoW05.fasta | samtools sort -@ 8 -O BAM -o glysoW05.to_minipan.bam
minipileup -f minilinerpan.fa -v -c -C -q 0 -s 0 -a 0 -l 1000 Wm82.to_minipan.bam > Wm82.to_minipan.PAV.vcf
minipileup -f minilinerpan.fa -v -c -C -q 0 -s 0 -a 0 -l 1000 glysoW05.to_minipan.bam > glysoW05.to_minipan.PAV.vcf
bcftools norm -m -any -f minilinerpan.fa -Oz -o Wm82.to_minipan.PAV.norm.vcf.gz Wm82.to_minipan.PAV.vcf
bcftools norm -m -any -f minilinerpan.fa -Oz -o glysoW05.to_minipan.PAV.norm.vcf.gz glysoW05.to_minipan.PAV.vcf
```

Extract deletions at least 50 bp long:

```bash
bcftools view -H Wm82.to_minipan.PAV.norm.vcf.gz | awk 'BEGIN{FS=OFS="\t"} $5!="." && $5!="*" && $5!~/^</ && length($4)>length($5) {start=$2-1+length($5); end=$2-1+length($4); if(end-start>=50) print $1,start,end,$1":"start"-"end,end-start}' | sort -k1,1 -k2,2n > W05_ins.bed
bcftools view -H glysoW05.to_minipan.PAV.norm.vcf.gz | awk 'BEGIN{FS=OFS="\t"} $5!="." && $5!="*" && $5!~/^</ && length($4)>length($5) {start=$2-1+length($5); end=$2-1+length($4); if(end-start>=50) print $1,start,end,$1":"start"-"end,end-start}' | sort -k1,1 -k2,2n > W82_ins.bed
```

Retain genome-specific regions overlapping TE annotations:

```bash
bedtools intersect -u -a W05_ins.bed -b minilinerpan.TE.bed > W05_TE_variants.bed
bedtools intersect -u -a W82_ins.bed -b minilinerpan.TE.bed > W82_TE_variants.bed
```

---

## 4. ATAC-seq and ChIP-seq

### 4.1 Read mapping

Check read quality with FastQC and MultiQC. Use fastp only when trimming is needed.

```bash
fastqc -t 16 *.fq.gz
multiqc .
bowtie2-build --threads 24 minilinerpan.fa minipan
```

Use local alignment for ATAC-seq:

```bash
while read sample; do bowtie2 -x minipan -1 ${sample}_1.fq.gz -2 ${sample}_2.fq.gz -X 2000 --dovetail --no-mixed --no-discordant --very-sensitive-local -p 24 --no-unal | samtools view -b -q 30 -F 1804 | samtools sort -@ 8 -O BAM -o ${sample}.sort.bam; sambamba markdup --remove-duplicates --nthreads 8 ${sample}.sort.bam ${sample}.rmDup.bam; samtools index ${sample}.rmDup.bam; done < ATAC.list
```

Use end-to-end alignment for ChIP-seq:

```bash
while read sample; do bowtie2 -x minipan -1 ${sample}_1.fq.gz -2 ${sample}_2.fq.gz -X 2000 --dovetail --no-mixed --no-discordant --very-sensitive -p 24 --no-unal | samtools view -b -q 30 -F 1804 | samtools sort -@ 8 -O BAM -o ${sample}.sort.bam; sambamba markdup --remove-duplicates --nthreads 8 ${sample}.sort.bam ${sample}.rmDup.bam; samtools index ${sample}.rmDup.bam; done < ChIP.list
```

### 4.2 Peak calling

```bash
GENOME_SIZE=$(awk '{s+=$2} END{print s}' minilinerpan.fa.fai)
while read sample; do macs3 callpeak -g ${GENOME_SIZE} -t ${sample}.rmDup.bam -n ${sample} -f BAMPE --keep-dup all -q 0.05 --broad; done < ATAC.list
```

For ATAC-seq and broad histone marks, reproducible peaks can be defined by reciprocal overlap between replicates:

```bash
awk '$9>=1.30103' W05_ATAC_rep1_peaks.broadPeak > W05_ATAC_rep1.q05.broadPeak
awk '$9>=1.30103' W05_ATAC_rep2_peaks.broadPeak > W05_ATAC_rep2.q05.broadPeak
bedtools intersect -f 0.5 -r -u -a W05_ATAC_rep1.q05.broadPeak -b W05_ATAC_rep2.q05.broadPeak > W05_ATAC_reproducible.bed
```

Example for a narrow mark such as H3K27ac:

```bash
macs3 callpeak -g ${GENOME_SIZE} -t W05_H3K27ac_rep1.rmDup.bam -c W05_input_rep1.rmDup.bam W05_input_rep2.rmDup.bam -n W05_H3K27ac_rep1 -f BAMPE --keep-dup all -q 0.05
idr --samples W05_H3K27ac_rep1_peaks.narrowPeak W05_H3K27ac_rep2_peaks.narrowPeak --input-file-type narrowPeak --rank p.value --idr-threshold 0.05 --output-file W05_H3K27ac_IDR.bed --plot
```

Example for a broad mark such as H3K27me3:

```bash
macs3 callpeak -g ${GENOME_SIZE} -t W05_H3K27me3_rep1.rmDup.bam -c W05_input_rep1.rmDup.bam W05_input_rep2.rmDup.bam -n W05_H3K27me3_rep1 -f BAMPE --keep-dup all -q 0.05 --broad
```

### 4.3 Peak quantification and differential analysis

Merge reproducible W05 and W82 peaks, then count reads in each sample with featureCounts:

```bash
cat W05_ATAC_reproducible.bed W82_ATAC_reproducible.bed | sort -k1,1 -k2,2n | bedtools merge -i - > All_ATACpeaks.bed
awk 'BEGIN{OFS="\t"; print "GeneID","Chr","Start","End","Strand"} {print $1":"$2+1"-"$3,$1,$2+1,$3,"."}' All_ATACpeaks.bed > All_ATACpeaks.saf
featureCounts -T 16 -p --countReadPairs -F SAF -a All_ATACpeaks.saf --fracOverlap 0.2 -o ATACpeaks_counts.tsv W05_ATAC_rep1.rmDup.bam W05_ATAC_rep2.rmDup.bam W82_ATAC_rep1.rmDup.bam W82_ATAC_rep2.rmDup.bam
```

Build a common peak set for each histone mark in the same way. Differential analysis should use raw counts, with `padj <= 0.05` as a possible cutoff.

---

## 5. BS-seq

```bash
mkdir methyindex
ln -s ../minilinerpan.fa methyindex/minilinerpan.fa
bismark_genome_preparation --bowtie2 --parallel 16 methyindex
while read sample; do bismark -N 0 -L 20 --bam --bowtie2 -p 8 --genome methyindex -1 ${sample}_1.fq.gz -2 ${sample}_2.fq.gz; deduplicate_bismark -bam -p ${sample}_1_bismark_bt2_pe.bam; bismark_methylation_extractor -p --parallel 4 --gzip --bedGraph --no_overlap --comprehensive --counts --CX_context --cytosine_report --genome_folder methyindex ${sample}_1_bismark_bt2_pe.deduplicated.bam; done < BS.list
```

Summarize CG, CHG, and CHH methylation separately. Convert the Bismark coverage file to BED and calculate methylated and unmethylated read counts within each TE variant region:

```bash
zcat W05_BS_rep1.bismark.cov.gz | awk 'BEGIN{OFS="\t"} {print $1,$2-1,$3,$5,$6}' > W05_BS_rep1.mCuC.bed
bedtools map -a W05_TE_variants.bed -b W05_BS_rep1.mCuC.bed -c 4,5 -o sum,sum > W05_BS_rep1.W05_TE_counts.bed
```

Process all samples and replicates in the same way, then compare regional methylation levels in R.

---

## 6. RNA-seq

```bash
gffread minilinerpan.gff3 -T -o minilinerpan.gtf
mkdir minipan_STAR_index
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir minipan_STAR_index --genomeFastaFiles minilinerpan.fa --sjdbGTFfile minilinerpan.gtf --sjdbOverhang 149 --genomeSAindexNbases 13
while read sample; do STAR --genomeDir minipan_STAR_index --runThreadN 16 --readFilesIn ${sample}_1.fq.gz ${sample}_2.fq.gz --readFilesCommand zcat --outFileNamePrefix ${sample}. --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --quantMode GeneCounts --outFilterScoreMinOverLread 0.5 --outFilterMatchNminOverLread 0.5 --outFilterIntronMotifs RemoveNoncanonical; done < RNAseq.list
```

Merge the `ReadsPerGene.out.tab` files and perform differential expression analysis with DESeq2 or TBtools. Use raw counts for differential analysis and TPM mainly for visualization.

---

## 7. Multi-omics integration

Intersect TE variants with differential ATAC-seq, histone modification, and methylation regions:

```bash
bedtools intersect -u -a W05_TE_variants.bed -b ATAC_DEpeaks.bed > W05_TE_ATAC.bed
bedtools intersect -u -a W05_TE_variants.bed -b H3K27ac_DEpeaks.bed > W05_TE_H3K27ac.bed
bedtools intersect -u -a W05_TE_variants.bed -b H3K27me3_DEpeaks.bed > W05_TE_H3K27me3.bed
```

Process W82-specific TE variants in the same way. Differentially expressed genes and TE-eQTLs can then be added using genomic overlap or linked gene IDs.

Regions derived from Wm82 can be converted from mini-pangenome coordinates back to Wm82 coordinates:

```bash
chainSwap Wm82_to_minipan.chain minipan_to_Wm82.chain
liftOver All_variation_epi.bed minipan_to_Wm82.chain All_variation_epi.on_Wm82.bed All_variation_epi.unmapped.bed
bedtools intersect -wa -wb -a All_variation_epi.on_Wm82.bed -b All_TAPTIP.bed > variation_epi_to_TAPTIP.bed
bedtools intersect -wa -wb -a variation_epi_to_TAPTIP.bed -b TE_eQTL.bed > variation_epi_to_TEeQTL.bed
```

W05-specific insertions do not have corresponding coordinates in Wm82. These regions should remain in mini-pangenome coordinates and be matched directly to `W05_ins.bed`.

---
