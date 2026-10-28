# CLL-CMA

Analysis code for the CLL Circulating Microenvironment Atlas (CLL-CMA): a disease-specific single-cell reference, bulk deconvolution of the CLL-map cohort, and associations of inferred cell fractions with survival, somatic alterations, splicing outliers, and germline variants.

This repository contains scripts only. Count matrices, clinical tables, splicing results, and genotypes are not included.

## How to run

Scripts look for inputs under `CLLCMA_alpha/` relative to the project root. Set that root before running them:

```bash
export CLL_TME_DIR=/path/to/CLLCMA_repo
cd "$CLL_TME_DIR"
```

R scripts that do not read `CLL_TME_DIR` fall back to `~/research/CLL_TME`. Python association scripts resolve the root from `CLL_TME_DIR`, or from the script location (`parents[2]` of `CLLCMA_alpha/scripts/`), which is this repository.

Place downloaded data in the paths named below. Do not commit those files.

## Pipeline

### 1. Single-cell atlas

| Step | Script | Output |
| --- | --- | --- |
| Harmonize cohorts | `CLLCMA_alpha/scripts/step_1_harmonize_datasets.R` | `CLLCMA_alpha/integrated_raw/` |
| Per-sample QC | `CLLCMA_alpha/scripts/step_1.5_qc_per_sample.R` | `integrated_raw/qc/` |
| Harmony integration and AIFI label transfer | `CLLCMA_alpha/scripts/step_2_annotate_integrated.R` | `integrated_annotated/` |
| Malignant B-cell assignment | `CLLCMA_alpha/scripts/step_3_assign_malignant_by_patient.R` | `seu_integrated_final.rds` |
| Deconvolution reference | `CLLCMA_alpha/scripts/step_4_export_bayesprism_ref.R` | `integrated_annotated/bayesprism/` |

Cohorts read by step 1, from `CLLCMA_alpha/datasets/<cohort>/`: PBMCpedia_HC, GSE111014 (day 0), GSE132065, GSE152469, GSE161610, GSE205704, GSE261917, GSE295489, GSE295491, Nadeu_2022. Each cohort needs a counts RDS and a cell-metadata TSV. Column maps are in `step_1_harmonize_datasets.R`.

Atlas figures: `plot_umap_predicted_cell_type.R`, `plot_cell_type_pies.R`, `plot_transfer_summary.R`.

### 2. B-cell stage programs

| Script | Role |
| --- | --- |
| `step_5_bcell_stage_deseq2.R` | Donor-stage pseudobulk DESeq2 (normal B, MBL, CLL, RT) |
| `plot_bcell_stage_markers.R`, `plot_bcell_stage_transition_volcanoes.R`, `replot_bcell_stage_heatmap_no_igv.R` | Stage figures |
| `gsea_bcell_stage_transitions.py` | Preranked GO / Reactome GSEA |
| `run_bcell_stage_infercnv.R` | inferCNV on B-lineage cells |
| `summarize_infercnv_by_chromosome.R`, `plot_infercnv_stage_ordered_heatmap.R`, `plot_infercnv_chromosome_summary.py`, `plot_infercnv_chromosome_arm_summary.py` | Chromosome-level summaries |

`summarize_infercnv_by_chromosome.R` expects `CLLCMA_alpha/references/hg38_chromosome_arm_boundaries.tsv`.

### 3. Bulk deconvolution and survival

| Script | Role |
| --- | --- |
| `deconv_cll_broad_2022.R` | quanTIseq, xCell, MCP-counter, CIBERSORT, and BayesPrism on CLL-map TPM |
| `cibersort_heatmap.R` | Fraction heatmaps (also used for BayesPrism) |
| `deconv_cox_forest.R` | Cox models of cell fraction versus OS and FFS |
| `cibersort_median_km_survival.R` | Median-split Kaplan–Meier curves |
| `performance_test.R` | 10-fold C-index and time-dependent AUC |
| `deepsurv_cv.py` | DeepSurv learner called by `performance_test.R` |

Bulk inputs:

- `CLLCMA_alpha/bulk/cll_broad_2022/data_mrna_seq_tpm.txt`
- `CLLCMA_alpha/bulk/cll_broad_2022/data_clinical_sample.txt`
- `CLLCMA_alpha/bulk/cll_broad_2022/data_clinical_patient.txt`

### 4. Somatic alterations and the microenvironment

| Script | Role |
| --- | --- |
| `genomic_microenv_association.R` | Mutations, arm CNVs, and focal CNVs versus BayesPrism fractions |
| `mutation_microenv_association.R` | Mutation-focused association and IGHV-adjusted models |
| `plot_alteration_microenv_heatmap.R` | Cliff’s delta heatmap |
| `alteration_tme_km_survival.R` | Survival splits for alteration-linked fractions |
| `rerun_associations_all_methods.sh` | Rerun mutation and genomic association scripts |

Also expects `data_mutations.txt`, `data_focal_cna.txt`, `data_armlevel_cna.txt`, and `data_sv.txt` in `CLLCMA_alpha/bulk/cll_broad_2022/`.

### 5. Splicing outliers and inferred fractions

| Script | Role |
| --- | --- |
| `splicing_outlier_tme_association.py` | Event outliers, gene 0/1 collapse, burden and per-gene models |
| `plot_splicing_outlier_tme.R` | Association figures |
| `splicing_outlier_rank_enrichment.py` | Rank by Cliff’s δ × −log10 adjusted P; ORA |
| `plot_rank_gsea_combined.py` | Rank-plus-enrichment figures |
| `splicing_outlier_cellmarker_gsea.py` | CellMarker preranked GSEA |
| `splicing_outlier_lr_prioritization.py` | Ligand–receptor prioritization |
| `plot_splicing_outlier_lr_networks.py` | LR network figures |

Splicing inputs live under `CLLCMA_alpha/bulk/cll_broad_2022_splicing/`. CellMarker and OmniPath/CellPhoneDB files are read from `CLLCMA_alpha/references/`. Set `DECONV_METHOD=bayesprism` (the default).

### 6. Germline association (phs000879)

Run from `CLLCMA_alpha/bulk/phs000879_GWAS/` after placing the approved PLINK genotype files there (`phs000879.bed/.bim/.fam`). Do not commit those files.

| Script | Role |
| --- | --- |
| `prepare_bayesprism_gwas_inputs.py` | Blom inverse-normal phenotypes and sex, age/10, IGHV covariates |
| `run_gemma_bayesprism_gwas.sh` | GEMMA relatedness matrix and LMM for 22 BayesPrism fractions |
| `plot_single_manhattan.R` | Manhattan and QQ plots (`topr`) |

`run_gemma_cibersort_gwas.sh` and `run_plink_cibersort_gwas.sh` are earlier CIBERSORT runs. The manuscript uses the BayesPrism GEMMA script.

External binaries: PLINK and GEMMA 0.98.3. The shell script activates a conda environment named `gwas`.

## Software

R packages used across the pipeline include Matrix, Seurat, harmony, data.table, ggplot2, survival, glmnet, ranger, DESeq2, immunedeconv, BayesPrism, infercnv, readxl, and topr.

Python packages are listed in `requirements-python.txt`. `performance_test.R` calls `deepsurv_cv.py` with PyTorch.

## What not to upload

- dbGaP genotypes, association scans, and subject-level clinical tables
- Single-cell count matrices and Seurat objects
- Bulk TPM, mutation, and splicing matrices
