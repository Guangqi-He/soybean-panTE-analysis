# Pan-core TE family analysis

This document describes the workflow used to estimate pan/core TE family curves and classify TE families into core, softcore, dispensable, and private categories.

## Purpose

This analysis summarizes the distribution of full-length TE families across soybean genomes and evaluates whether each TE family is shared by all, most, several, or only one accession.

## Required tools

- EDTA utility script `find_flTE.pl`
- R with the `tidyverse` package
- awk, grep, sort

## Required input files

- `Wm82.fasta.mod.out`: RepeatMasker `.out` file generated from pan-genome TE library annotation
- `all_full_length_TE_families.tsv`: a two-column table containing sample ID and full-length TE family ID
- `all_TE_family_presence.tsv`: a two-column table containing sample ID and TE family ID
- `all_TE_positions.bed`: TE annotation table with at least five columns: chromosome, start, end, TE family ID, and TE type
- `sample_metadata.tsv`: sample metadata table with a `sample` column

## Step 1. Extract full-length TE records from one genome

```bash
perl find_flTE.pl Wm82.fasta.mod.out \
  | grep -v -P 'CL569186.1|AF013103.1|\)n|cent|Cent|telo|knob|TR-1|osed|sela|A-rich|G-rich' \
  > Wm82.full_length_TE.txt
```

## Step 2. Obtain a non-redundant full-length TE family list

In the EDTA full-length TE output, the TE family ID is typically stored in column 10.

```bash
awk '{print $10}' Wm82.full_length_TE.txt | sort -u > Wm82.full_length_TE_family.list

awk 'BEGIN{OFS="\t"} {print "Wm82", $1}' \
  Wm82.full_length_TE_family.list \
  > Wm82.full_length_TE_family.tsv
```

For the full analysis, equivalent two-column tables from all genomes were combined into `all_full_length_TE_families.tsv`.

## Step 3. Estimate pan/core TE family curves

The following R script uses a combined two-column table with columns `sample` and `family`. It uses `replicate()` and `sapply()` to avoid nested shell loops.

```r
library(tidyverse)

# Input table: sample family
te <- read_tsv(
  "all_full_length_TE_families.tsv",
  col_names = c("sample", "family"),
  show_col_types = FALSE
) %>%
  distinct(sample, family)

sample_list <- sort(unique(te$sample))
sample_count <- length(sample_list)
n_bootstrap <- 1000
set.seed(123)

one_bootstrap <- function() {
  pan_count <- integer(sample_count)
  core_count <- integer(sample_count)

  for (k in seq_len(sample_count)) {
    selected_samples <- sample(sample_list, size = k, replace = FALSE)
    selected_te <- filter(te, sample %in% selected_samples)

    family_freq <- table(selected_te$family)
    pan_count[k] <- length(family_freq)
    core_count[k] <- sum(family_freq == k)
  }

  list(pan = pan_count, core = core_count)
}

boot <- replicate(n_bootstrap, one_bootstrap(), simplify = FALSE)

pan_matrix <- do.call(rbind, lapply(boot, `[[`, "pan"))
core_matrix <- do.call(rbind, lapply(boot, `[[`, "core"))

write.table(pan_matrix, "panTE_bootstrap1000.tsv", sep = "\t", quote = FALSE)
write.table(core_matrix, "coreTE_bootstrap1000.tsv", sep = "\t", quote = FALSE)
```

## Step 4. Classify All TE families by their frequency across genomes

This example assumes 47 genomes. Adjust `n_genomes` if a different number of genomes is used.

```r
library(tidyverse)

n_genomes <- 47

presence <- read_tsv(
  "all_TE_family_presence.tsv",
  col_names = c("sample", "family"),
  show_col_types = FALSE
) %>%
  distinct(sample, family) %>%
  filter(!str_detect(family, fixed(")n"))) %>%
  filter(!str_detect(family, regex("rich", ignore_case = TRUE)))

family_class <- presence %>%
  count(family, name = "n_samples") %>%
  mutate(class = case_when(
    n_samples == n_genomes ~ "core",
    n_samples >= n_genomes - 2 & n_samples < n_genomes ~ "softcore",
    n_samples == 1 ~ "private",
    n_samples >= 2 & n_samples <= n_genomes - 3 ~ "dispensable",
    TRUE ~ NA_character_
  ))

write_tsv(family_class, "TE_family_frequency_class.tsv")
```

## Step 5. Summarize TE copy number and total length by family class

```r
library(tidyverse)

te_pos <- read_tsv(
  "all_TE_positions.bed",
  col_names = c("chr", "start", "end", "family", "type"),
  show_col_types = FALSE
) %>%
  mutate(length = end - start + 1)

family_class <- read_tsv("TE_family_frequency_class.tsv", show_col_types = FALSE)
sample_meta <- read_tsv("sample_metadata.tsv", show_col_types = FALSE)

te_summary <- te_pos %>%
  left_join(family_class, by = "family") %>%
  filter(!is.na(class)) %>%
  group_by(class) %>%
  summarise(
    copy_number = n(),
    total_length_Mb = sum(length) / 1e6,
    family_number = n_distinct(family),
    .groups = "drop"
  )

write_tsv(te_summary, "TE_family_class_summary.tsv")
```

## Notes

- `core` families are present in all genomes.
- `softcore` families are absent from one or two genomes.
- `dispensable` families are absent from more than two genomes but present in at least two genomes.
- `private` families are present in only one genome.
- Remove tandem, centromeric, telomeric, and low-complexity entries before pan/core curve estimation when these records are not considered true TE families.
