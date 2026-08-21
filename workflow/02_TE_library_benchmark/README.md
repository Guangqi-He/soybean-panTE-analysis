# CACTA Annotation Assessment on the Wm82.a1 Genome

This workflow evaluates the unusually high CACTA abundance obtained when the SoyTEdb CACTA library is used alone, tests the contribution of nested sequences, and examines how the complete SoyTEdb library reclassifies the excess CACTA-only matches. The final section benchmarks the pan-genome TE library against the complete SoyTEdb annotation.

## Input files

```text
Wm82.a1.fa
soyTEdb.fa
soyTEdb_CACTA.fa
Glycine_panTElib.fa
```

The following scripts are assumed to be available in the working environment:

```text
cleanup_tandem.pl
cleanup_nested.pl
count_base.pl
buildSummary.pl
lib-test.py
rmout2bed.py
classify_rm_overlap.py
compare_CACTA_SW.py
```

Set the main variables once before running the workflow:

```bash
GENOME=Wm82.a1.fa
FULL_LIB=soyTEdb.fa
RAW_CACTA_LIB=soyTEdb_CACTA.fa
CLEAN_CACTA_LIB=soyTEdb_CACTA_cleanNest.fa
PAN_LIB=Glycine_panTElib.fa

FULL_OUT=Wm82_a1.SoyTEdb.out
RAW_OUT=Wm82_a1.SoyTEdb_CACTA.out
CLEAN_OUT=Wm82_a1.SoyTEdb_CACTA_cleanNest.out
PAN_OUT=Wm82_a1.panTElib.out

FULL_BED=Wm82_a1.SoyTEdb.bed
RAW_BED=Wm82_a1.SoyTEdb_CACTA.bed
CLEAN_BED=Wm82_a1.SoyTEdb_CACTA_cleanNest.bed

FULL_CACTA=Wm82_a1.SoyTEdb.full_CACTA.merge.bed
RAW=Wm82_a1.SoyTEdb_CACTA.merge.bed
CLEAN=Wm82_a1.SoyTEdb_CACTA_cleanNest.merge.bed

FAI=${GENOME}.fai
THREADS=48
```

## 1. Reproduce the CACTA abundance difference

### 1.1 RepeatMasker annotation with the complete SoyTEdb library

```bash
RepeatMasker -e ncbi -pa ${THREADS} -q -no_is -norna -nolow -div 40 -lib ${FULL_LIB} -cutoff 225 ${GENOME}
mv ${GENOME}.out ${FULL_OUT}
```

### 1.2 RepeatMasker annotation with the CACTA-only SoyTEdb library

```bash
RepeatMasker -e ncbi -pa ${THREADS} -q -no_is -norna -nolow -div 40 -lib ${RAW_CACTA_LIB} -cutoff 225 ${GENOME}
mv ${GENOME}.out ${RAW_OUT}
```

Using the complete SoyTEdb library gives approximately 2% CACTA content, whereas using the isolated CACTA library gives approximately 10%.

## 2. Inspect nested sequences in SoyTEdb CACTA models

TE-Aid can be used to inspect representative CACTA sequences before cleanup.

```bash
seqkit split -i -O CACTA_split ${RAW_CACTA_LIB}

for fa in CACTA_split/*.fa
do
    TE-Aid -q "$fa" -g ${GENOME} -t -title before -v_x_line_1 0 -v_x_line_2 0 -o TE_Aid_before
done
```

This step is intended for visual inspection rather than genome-wide quantification.

## 3. Generate a nested-sequence-cleaned CACTA library

First remove tandem/repetitive artifacts from the complete SoyTEdb library and then apply the EDTA nested-sequence cleanup.

```bash
cleanup_tandem.pl -f ${FULL_LIB} -misschar N -nc 50000 -nr 0.9 -minlen 100 -minscore 3000 -trf 1 -cleanN 1 > soyTEdb_clean_tandem.fa
cleanup_nested.pl -in soyTEdb_clean_tandem.fa -threads ${THREADS}
seqkit grep -r -p "DNA/CACTA" soyTEdb_clean_tandem.fa.cln > ${CLEAN_CACTA_LIB}
```

## 4. Annotate the genome with the cleaned CACTA library

```bash
RepeatMasker -e ncbi -pa ${THREADS} -q -no_is -norna -nolow -div 40 -lib ${CLEAN_CACTA_LIB} -cutoff 225 ${GENOME}
mv ${GENOME}.out ${CLEAN_OUT}
```

The cleaned CACTA-only library typically reduces the apparent CACTA abundance from approximately 10% to approximately 5%.

## 5. Generate RepeatMasker summary tables

```bash
count_base.pl ${GENOME} > Wm82_a1.stats
buildSummary.pl -maxDiv 40 -stats Wm82_a1.stats ${FULL_OUT} > Wm82_a1.SoyTEdb.sum 2>/dev/null
buildSummary.pl -maxDiv 40 -stats Wm82_a1.stats ${RAW_OUT} > Wm82_a1.SoyTEdb_CACTA.sum 2>/dev/null
buildSummary.pl -maxDiv 40 -stats Wm82_a1.stats ${CLEAN_OUT} > Wm82_a1.SoyTEdb_CACTA_cleanNest.sum 2>/dev/null
```

## 6. Convert RepeatMasker output to BED and merge intervals

```bash
python3 rmout2bed.py ${FULL_OUT} > ${FULL_BED}
python3 rmout2bed.py ${RAW_OUT} > ${RAW_BED}
python3 rmout2bed.py ${CLEAN_OUT} > ${CLEAN_BED}
```

Extract CACTA annotations from the complete SoyTEdb result:

```bash
awk '$5=="DNA/CACTA"' ${FULL_BED} | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i - > ${FULL_CACTA}
```

Merge CACTA-only annotations:

```bash
cut -f1-3 ${RAW_BED} | sort -k1,1 -k2,2n | bedtools merge -i - > ${RAW}
cut -f1-3 ${CLEAN_BED} | sort -k1,1 -k2,2n | bedtools merge -i - > ${CLEAN}
```

## 7. Recalculate genome coverage

Create the FASTA index if it is not already available:

```bash
samtools faidx ${GENOME}
```

Calculate genome size and the non-overlapping CACTA coverage for each annotation set:

```bash
GENOME_SIZE=$(awk '{sum+=$2} END{print sum}' ${FAI})
echo "Genome size: ${GENOME_SIZE}"

for f in ${FULL_CACTA} ${RAW} ${CLEAN}
do
    BP=$(awk '{sum+=$3-$2} END{print sum+0}' "$f")
    PCT=$(awk -v bp=${BP} -v genome=${GENOME_SIZE} \
      'BEGIN{printf "%.6f",bp/genome*100}')
    echo -e "$f\t${BP}\t${PCT}"
done
```

Example result:

```text
Wm82_a1.SoyTEdb.full_CACTA.merge.bed              12813022    1.316391
Wm82_a1.SoyTEdb_CACTA.merge.bed                  106908993   10.983676
Wm82_a1.SoyTEdb_CACTA_cleanNest.merge.bed         53056919    5.450991
```

## 8. Basic consistency checks

Check whether the cleaned CACTA annotation contains intervals absent from the raw CACTA annotation:

```bash
bedtools subtract -a ${CLEAN} -b ${RAW} > QC.clean_not_raw.bed
awk '{sum+=$3-$2} END{print sum+0}' QC.clean_not_raw.bed
```

Compare interval overlap:

```bash
bedtools jaccard -a ${RAW} -b ${CLEAN}
bedtools jaccard -a ${CLEAN} -b ${FULL_CACTA}
```

Observed values:

```text
RAW vs CLEAN:
intersection    union        jaccard    n_intersections
52993922        106971990    0.4954     150795

CLEAN vs FULL_CACTA:
intersection    union       jaccard     n_intersections
13839336        53275471    0.259769    26859
```

## 9. Define CACTA annotation difference sets

### 9.1 Cleanup-sensitive regions

These regions are detected as CACTA only when the uncleaned CACTA library is used and disappear after nested-sequence cleanup.

```bash
bedtools subtract -a ${RAW} -b ${CLEAN} | sort -k1,1 -k2,2n | bedtools merge -i - > CACTA.cleanup_sensitive.bed
awk '{sum+=$3-$2} END{print sum+0}' CACTA.cleanup_sensitive.bed
```

Observed total:

```text
53915071 bp
```

### 9.2 Residual CACTA-only regions after cleanup

These regions remain detectable with the cleaned CACTA-only library but are not annotated as CACTA when all SoyTEdb TE families compete during RepeatMasker annotation.

```bash
bedtools subtract -a ${CLEAN} -b ${FULL_CACTA} | sort -k1,1 -k2,2n | bedtools merge -i - > CACTA.cleaned_residual.bed
awk '{sum+=$3-$2} END{print sum+0}' CACTA.cleaned_residual.bed
```

Observed total:

```text
39217583 bp
```

### 9.3 Concordant CACTA regions

These regions are consistently detected by both the cleaned CACTA-only library and the complete SoyTEdb library.

```bash
bedtools intersect -a ${CLEAN} -b ${FULL_CACTA} | sort -k1,1 -k2,2n | bedtools merge -i - > CACTA.cleaned_concordant.bed
awk '{sum+=$3-$2} END{print sum+0}' CACTA.cleaned_concordant.bed
```

Observed total:

```text
13839336 bp
```

## 10. Determine how the complete SoyTEdb library annotates the difference sets

Intersect each difference set with the original, unmerged complete SoyTEdb annotation so that the RepeatMasker annotation fields and SW scores are retained.

```bash
bedtools intersect -a CACTA.cleanup_sensitive.bed -b ${FULL_BED} -wao > cleanup_sensitive_vs_fullSoyTEdb.tsv
bedtools intersect -a CACTA.cleaned_residual.bed -b ${FULL_BED} -wao > cleaned_residual_vs_fullSoyTEdb.tsv
```

Because multiple SoyTEdb annotations can overlap the same region, retain the overlapping annotation with the highest RepeatMasker SW score:

```bash
python3 classify_rm_overlap.py --input cleanup_sensitive_vs_fullSoyTEdb.tsv --genome-size ${GENOME_SIZE} --prefix cleanup_sensitive
python3 classify_rm_overlap.py --input cleaned_residual_vs_fullSoyTEdb.tsv --genome-size ${GENOME_SIZE} --prefix cleaned_residual
```

## 11. Compare RepeatMasker SW scores

For cleanup-sensitive regions, compare the best CACTA-only match from the raw CACTA library with the best annotation from the complete SoyTEdb library:

```bash
python3 compare_CACTA_SW.py --regions CACTA.cleanup_sensitive.bed --cacta-bed ${RAW_BED} --full-bed ${FULL_BED} --label Cleanup-sensitive --out CACTA.cleanup_sensitive_SW --min-overlap 80 --min-region-cov 0.5
```

For residual regions, use the cleaned CACTA annotation:

```bash
python3 compare_CACTA_SW.py --regions CACTA.cleaned_residual.bed --cacta-bed ${CLEAN_BED} --full-bed ${FULL_BED} --label Residual-after-cleanup --out CACTA.cleaned_residual_SW --min-overlap 80 --min-region-cov 0.5
```

Combine paired results:

```bash
head -n 1 CACTA.cleanup_sensitive_SW.paired.tsv > CACTA.excess_SW.all.tsv
tail -n +2 CACTA.cleanup_sensitive_SW.paired.tsv >> CACTA.excess_SW.all.tsv
tail -n +2 CACTA.cleaned_residual_SW.paired.tsv >> CACTA.excess_SW.all.tsv
```

### Plot SW-score comparisons in R

```r
library(ggplot2)

dat <- read.table("CACTA.excess_SW.all.tsv",header = TRUE,sep = "\t",check.names = FALSE)
dat <- dat[!is.na(dat$CACTA_SW) & !is.na(dat$Full_SW), ]

p <- ggplot(dat, aes(x = CACTA_SW, y = Full_SW, color = Full_broad_class)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_point(alpha = 0.55, size = 1.5) +
    facet_wrap(~ region_set) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
        x = "Best SW score with CACTA-only library",
        y = "Best non-CACTA SW score with complete SoyTEdb",
        color = "Best annotation in complete SoyTEdb"
    ) +
    theme_classic()

ggsave("CACTA.excess_SW.scatter.pdf",p,width = 7,height = 4.5)
with(dat, tapply(Delta_SW > 0, region_set, mean))
```
Observed proportions:

```text
Cleanup-sensitive        0.9934628
Residual-after-cleanup   0.9815651
```

Thus, 99.35% of cleanup-sensitive regions and 98.16% of residual CACTA-only regions had higher SW scores for non-CACTA TE models in the complete SoyTEdb library than for CACTA models in the corresponding CACTA-only annotation.

## 12. Benchmark the pan-genome TE library

Annotate Wm82.a1 with the complete pan-genome TE library:

```bash
RepeatMasker -e ncbi -pa ${THREADS} -q -no_is -norna -nolow -div 40 -lib ${PAN_LIB} -cutoff 225 ${GENOME}
mv ${GENOME}.out ${PAN_OUT}
```

Compare the pan-genome TE annotation with the complete SoyTEdb annotation using EDTA `lib-test.py`:

```bash
python3 lib-test.py --genome ${GENOME} --test ${PAN_OUT} --reference ${FULL_OUT} --extended_report
```

This command evaluates the complete pan-genome TE annotation against the complete SoyTEdb annotation. If a CACTA-only benchmark is required, CACTA-only RepeatMasker outputs should be generated for both libraries before running `lib-test.py`.
