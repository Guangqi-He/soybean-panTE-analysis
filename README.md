# Soybean pan-genome transposable element analysis

This repository contains relevant data, custom scripts, and workflow notes for the manuscript:

**A pan-genome transposable element library reveals TE-derived variants and epigenetic rewiring of gene expression during soybean domestication**

Transposable elements (TEs) are key drivers of genome evolution, but their contribution to soybean domestication has remained incompletely understood. In this study, we constructed a comprehensive *Glycine* pan-genome TE library from 53 high-quality genomes and identified 18,797 non-redundant TE families. We further investigated TE presence/absence polymorphisms, TE-derived variants, TE-associated epigenetic changes, and their effects on gene expression during soybean domestication.

This repository provides selected example data, custom scripts, and step-by-step descriptions for the major analyses performed in the study.

## Repository structure
- [01_pan_genome_TE_library_construction](https://github.com/Guangqi-He/soybean-panTE-analysis/tree/main/workflow/01_pan_genome_TE_library_construction): Construction of pan-genome TE library.
- [02_TE_library_benchmark](https://github.com/Guangqi-He/soybean-panTE-analysis/tree/main/workflow/02_TE_library_benchmark): Comparison between different TE library.
- [03_pan_core_TE_family_analysis](https://github.com/Guangqi-He/soybean-panTE-analysis/tree/main/workflow/03_pan_core_TE_family_analysis): Classification of TE families.
- [04_RT_domain_phylogenetic_analysis](https://github.com/Guangqi-He/soybean-panTE-analysis/tree/main/workflow/04_RT_domain_phylogenetic_analysis): phylogenetic reconstruction for selected retrotransposon families.
- [05_TE_distance_to_TSS_and_genes](https://github.com/Guangqi-He/soybean-panTE-analysis/tree/main/workflow/05_TE_distance_to_TSS_and_genes): Characterization of TE distributions relative to genes and TSSs.
- [06_syntenic_diversity](https://github.com/Guangqi-He/soybean-panTE-analysis/tree/main/workflow/06_syntenic_diversity): Quantification of syntenic diversity.
- [07_TAP_TIP_identification](https://github.com/Guangqi-He/soybean-panTE-analysis/tree/main/workflow/07_TAP_TIP_identification): Identification of TE absence and insertion polymorphisms.
- [08_eQTL_analysis](https://github.com/Guangqi-He/soybean-panTE-analysis/tree/main/workflow/08_eQTL_analysis): eQTL analysis of TAPs/TIPs.
- [09_epigenetic_analysis](https://github.com/Guangqi-He/soybean-panTE-analysis/tree/main/workflow/09_epigenetic_analysis): Integration of TE-derived variants with epigenetic data.
- [10_TE_expression_analysis](https://github.com/Guangqi-He/soybean-panTE-analysis/tree/main/workflow/10_TE_expression_analysis): Quantification of TE family expression.
- `data/`: Selected example datasets and processed tables used to demonstrate the analysis workflow.
- `scripts/`: Custom scripts used for TE-related analyses.

## Usage notes

Most scripts were designed for command-line execution in a Linux environment. Users should check file paths, software versions, and input formats before running the scripts on their own datasets.

Example usage and required input files are provided in the corresponding subdirectories when applicable.

## Contact

For questions about this repository, please contact GuangKi (guangqi@zju.edu.cn).
