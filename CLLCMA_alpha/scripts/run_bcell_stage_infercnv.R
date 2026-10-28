#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(biomaRt)
  library(infercnv)
})

project_dir <- Sys.getenv("CLL_TME_DIR", unset = path.expand("~/research/CLL_TME"))
setwd(normalizePath(project_dir, mustWork = TRUE))

seu_file <- "CLLCMA_alpha/integrated_annotated/seu_integrated_final.rds"
base_res_dir <- "CLLCMA_alpha/results/bcell_stage_deseq2"
out_dir <- file.path(base_res_dir, "publication_figures", "infercnv_bcell_stage")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(42L)
max_cells_per_stage <- as.integer(Sys.getenv("INFERCNV_MAX_CELLS_PER_STAGE", "3000"))
max_cells_per_donor_stage <- as.integer(Sys.getenv("INFERCNV_MAX_CELLS_PER_DONOR_STAGE", "250"))
num_threads <- as.integer(Sys.getenv("INFERCNV_THREADS", "8"))

stage_levels <- c("HC_B", "MBL", "CLL", "RT")
stage_labels <- c(HC_B = "Normal_B", MBL = "MBL", CLL = "CLL", RT = "RT")
hc_b_labels <- c("Naive B cell", "Memory B cell", "Transitional B cell", "Effector B cell")
nadeu_rs_like <- c(
  "RS", "RT", "RT proliferative", "RT quiescent",
  "TP53INP1hi RT", "MIR155HGhi RT", "CCND2hi RT", "CCND2lo RT",
  "CXCR4loCD27hi RT"
)
nadeu_cll_like <- c(
  "CLL", "CXCR4hiCD27lo", "CXCR4loCD27hi",
  "MIR155HGhi", "MIR155HGhi CLL",
  "CD83hiMIR155HGhi", "CD83loMIR155HGhi",
  "MZB1hiIGHMhiXBP1hi"
)

assign_b_stage <- function(cell_type_final, disease, dataset) {
  ct <- as.character(cell_type_final)
  dis <- as.character(disease)
  ds <- as.character(dataset)
  stage <- rep(NA_character_, length(ct))
  stage[ct %in% hc_b_labels & dis == "HC"] <- "HC_B"
  stage[ct == "MBL"] <- "MBL"
  stage[ct == "CLL"] <- "CLL"
  stage[ct == "RT"] <- "RT"
  stage[ct %in% nadeu_cll_like & dis == "CLL"] <- "CLL"
  stage[ct %in% nadeu_rs_like] <- "RT"
  stage[dis == "RS" & ct %in% nadeu_cll_like] <- "RT"
  stage[ds == "PBMCpedia_HC" & !ct %in% hc_b_labels] <- NA_character_
  stage
}

get_counts <- function(seu) {
  tryCatch(
    LayerData(seu, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(seu, assay = "RNA", slot = "counts")
  )
}

balanced_sample <- function(dt, max_total, max_per_donor) {
  sampled <- dt[, {
    n <- min(.N, max_per_donor)
    .SD[sample.int(.N, n)]
  }, by = donor]
  if (nrow(sampled) > max_total) {
    sampled <- sampled[sample.int(nrow(sampled), max_total)]
  }
  sampled
}

message("Loading Seurat object: ", seu_file)
seu <- readRDS(seu_file)
meta <- as.data.table(seu[[]], keep.rownames = "cell_id")
need <- c("dataset", "sample", "donor", "disease", "cell_type_final")
miss <- setdiff(need, names(meta))
if (length(miss)) stop("Missing metadata columns: ", paste(miss, collapse = ", "))

meta[, b_stage := assign_b_stage(cell_type_final, disease, dataset)]
meta <- meta[!is.na(b_stage) & b_stage %in% stage_levels]
meta[, b_stage := factor(b_stage, levels = stage_levels)]
meta[, infercnv_group := stage_labels[as.character(b_stage)]]

message("Eligible B-lineage cells by stage before inferCNV sampling:")
print(meta[, .N, by = b_stage][order(b_stage)])

keep_meta <- rbindlist(lapply(stage_levels, function(st) {
  dt <- meta[b_stage == st]
  if (!nrow(dt)) stop("No cells available for stage ", st)
  balanced_sample(dt, max_cells_per_stage, max_cells_per_donor_stage)
}))
setorder(keep_meta, b_stage, dataset, donor, sample, cell_id)
message("Cells selected for inferCNV:")
print(keep_meta[, .(n_cells = .N, n_donors = uniqueN(donor), n_datasets = uniqueN(dataset)), by = .(b_stage, infercnv_group)])
fwrite(keep_meta, file.path(out_dir, "infercnv_selected_cells_metadata.tsv"), sep = "\t")

message("Extracting raw counts for selected cells ...")
counts <- get_counts(seu)
missing_cells <- setdiff(keep_meta$cell_id, colnames(counts))
if (length(missing_cells)) stop("Selected cells missing from counts: ", length(missing_cells))
counts <- counts[, keep_meta$cell_id, drop = FALSE]
rm(seu)
gc(verbose = FALSE)

coord_cache <- file.path(base_res_dir, "publication_figures", "cnv_heatmap", "gene_coordinates_biomart.tsv")
if (file.exists(coord_cache)) {
  coords <- fread(coord_cache)
} else {
  message("Retrieving gene coordinates from Ensembl BioMart ...")
  mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast")
  coords <- getBM(
    attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
    filters = "hgnc_symbol",
    values = rownames(counts),
    mart = mart
  )
  coords <- as.data.table(coords)
  setnames(coords, c("hgnc_symbol", "chromosome_name", "start_position", "end_position"), c("gene", "chr", "start", "end"))
}
coords <- coords[gene %in% rownames(counts)]
coords <- coords[gene != "" & as.character(chr) %in% as.character(1:22) & !is.na(start) & !is.na(end)]
coords[, chr := paste0("chr", as.character(chr))]
coords[, width := end - start + 1]
setorder(coords, gene, chr, start, -width)
coords <- coords[!duplicated(gene)]
coords[, chr_num := as.integer(sub("chr", "", chr))]
setorder(coords, chr_num, start, end, gene)
message("Genes retained with autosomal genomic coordinates: ", nrow(coords), " / ", nrow(counts))

counts <- counts[coords$gene, , drop = FALSE]
annotations <- keep_meta[, .(cell_id, infercnv_group)]
gene_order <- coords[, .(gene, chr, start, end)]

counts_file <- file.path(out_dir, "infercnv_raw_counts.rds")
anno_file <- file.path(out_dir, "infercnv_cell_annotations.tsv")
gene_file <- file.path(out_dir, "infercnv_gene_order.tsv")
saveRDS(counts, counts_file, compress = FALSE)
fwrite(annotations, anno_file, sep = "\t", col.names = FALSE)
fwrite(gene_order, gene_file, sep = "\t", col.names = FALSE)

message("Creating inferCNV object ...")
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts,
  annotations_file = anno_file,
  gene_order_file = gene_file,
  ref_group_names = "Normal_B",
  delim = "\t",
  min_max_counts_per_cell = c(100, Inf),
  chr_exclude = c("chrX", "chrY", "chrM")
)

run_dir <- file.path(out_dir, "infercnv_run")
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
message("Running inferCNV -> ", run_dir)
infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = run_dir,
  cluster_by_groups = TRUE,
  cluster_references = FALSE,
  denoise = TRUE,
  HMM = FALSE,
  num_threads = num_threads,
  output_format = "png",
  plot_steps = FALSE,
  save_rds = TRUE,
  save_final_rds = TRUE,
  no_prelim_plot = TRUE,
  no_plot = FALSE
)

saveRDS(infercnv_obj, file.path(out_dir, "infercnv_obj_final.rds"), compress = FALSE)

caption <- c(
  "Figure X. inferCNV analysis of B-cell disease stages using normal B cells as reference.",
  "Raw single-cell counts from CLL-CMA B-lineage populations were ordered by genomic coordinate and analyzed with inferCNV, using normal B cells as the reference group and MBL, CLL, and RT cells as observation groups.",
  paste0("Cells were subsampled in a donor-balanced manner with up to ", max_cells_per_stage, " cells per stage and up to ", max_cells_per_donor_stage, " cells per donor-stage stratum."),
  "Red and blue segments indicate inferred broad expression gains and losses, respectively, relative to the normal B-cell reference. inferCNV results represent RNA-expression-based copy-number inference and should be interpreted alongside matched genomic CNA data when available."
)
writeLines(caption, file.path(out_dir, "infercnv_bcell_stage_caption.txt"))
message("Done inferCNV output: ", normalizePath(out_dir))
