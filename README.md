# Soybean pan-genome transposable element analysis

This repository contains relevant data, custom scripts, and workflow notes for the manuscript:

**A pan-genome transposable element library reveals TE-derived variants and epigenetic rewiring of gene expression during soybean domestication**

Transposable elements (TEs) are key drivers of genome evolution, but their contribution to soybean domestication has remained incompletely understood. In this study, we constructed a comprehensive *Glycine* pan-genome TE library from 53 high-quality genomes and identified 18,797 non-redundant TE families. We further investigated TE presence/absence polymorphisms, TE-derived variants, TE-associated epigenetic changes, and their effects on gene expression during soybean domestication.

This repository provides selected example data, custom scripts, and step-by-step descriptions for the major analyses performed in the study.

## Repository structure

- `data/`  
  Selected example datasets and processed tables used to demonstrate the analysis workflow.

- `scripts/`  
  Custom scripts used for TE-related analyses, data processing, statistical analysis, and figure generation.

- `01_pan_genome_TE_library_construction/`  
  Workflow for constructing the *Glycine* pan-genome TE library, including TE library cleaning, genome-level TE annotation, pan-genome TE library generation, and reclassification of unknown TE families.

- `02_TE_library_benchmark/`  
  Scripts and notes for benchmarking the pan-genome TE library against existing soybean TE resources and evaluating annotation performance.

- `03_pan_core_TE_family_analysis/`  
  Analysis of core, softcore, dispensable, and private TE families across soybean accessions.

- `04_RT_domain_phylogenetic_analysis/`  
  Workflow for extracting reverse transcriptase (RT) domains and performing phylogenetic analysis of LTR retrotransposon families.

- `05_TE_distance_to_TSS_and_genes/`  
  Analysis of TE distribution relative to genes and transcription start sites (TSSs), including distance calculation and gene-proximal TE classification.

- `06_syntenic_diversity/`  
  Analysis of syntenic diversity and its relationship with TE density and gene density across the soybean genome.

- `07_TAP_TIP_identification/`  
  Identification and classification of TE absence polymorphisms (TAPs) and TE insertion polymorphisms (TIPs) from population-scale TE-derived variants.

- `08_epigenetic_analysis/`  
  Analysis of chromatin accessibility, DNA methylation, and TE-associated epigenetic rewiring.

- `09_TE_expression_analysis/`  
  TE expression quantification and analysis of transcriptionally expressed TE families across soybean tissues and accessions.

## Workflow overview

The analyses are organized in approximately the following order:

1. **Pan-genome TE library construction**  
   Construction of a non-redundant *Glycine* pan-genome TE library from multiple high-quality genomes.

2. **TE library benchmark**  
   Comparison of the pan-genome TE library with existing soybean TE resources and assessment of annotation consistency.

3. **Pan-core TE family analysis**  
   Classification of TE families into core, softcore, dispensable, and private categories.

4. **RT domain phylogenetic analysis**  
   Extraction of RT domains and phylogenetic reconstruction for selected retrotransposon families.

5. **TE location relative to genes**  
   Characterization of TE distributions relative to genes, TSSs, and regulatory regions.

6. **Syntenic diversity analysis**  
   Quantification of syntenic diversity and evaluation of its association with TE density and gene density.

7. **TAP/TIP identification**  
   Identification of TE absence and insertion polymorphisms across soybean populations.

8. **Epigenetic analysis**  
   Integration of TE-derived variants with chromatin accessibility and DNA methylation data.

9. **TE expression analysis**  
   Quantification of TE expression and identification of transcriptionally expressed TE families.

## Data availability

Large raw datasets, genome assemblies, and full-scale intermediate files are not stored directly in this repository because of file-size limitations. This repository mainly contains selected example data, processed tables, custom scripts, and workflow documentation.

The complete data sources and accession information are described in the manuscript and supplementary materials.

## Usage notes

Most scripts were designed for command-line execution in a Linux environment. Users should check file paths, software versions, and input formats before running the scripts on their own datasets.

Example usage and required input files are provided in the corresponding subdirectories when applicable.


## Contact

For questions about this repository, please contact the corresponding authors of the manuscript.
