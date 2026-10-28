#!/usr/bin/env Rscript

# Immune deconvolution of the Broad CLL bulk RNA-seq cohort using the
# CLLCMA_alpha scRNA BayesPrism / CIBERSORT reference.
#
# Inputs:
#   CLLCMA_alpha/bulk/cll_broad_2022/data_mrna_seq_tpm.txt
#   CLLCMA_alpha/integrated_annotated/bayesprism/
#
# Outputs:
#   CLLCMA_alpha/results/deconvolution/cll_broad_2022/
#     <method>_sample_by_feature.tsv/.rds
#     all_methods_long.tsv
#
# Methods: quanTIseq, xCell, MCP-counter,
# CIBERSORT (signature from CLLCMA GSE scRNA cell types),
# and BayesPrism (CLLCMA GSE scRNA counts; key = "CLL").
#
# BayesPrism allows a single malignant key ("CLL"). If an older reference still
# contains RS as a cell_type, those cells are collapsed into CLL while retaining
# RS_* cell_state labels. The current fine-type CLLCMA export uses CLL-only
# malignant cells (no RS/MBL in the malignant compartment).
#
# ESTIMATE_TYPE distinguishes compositional estimates ("fraction") from
# enrichment/abundance scores, which must not be interpreted as percentages.
#
# Bulk input is TPM; BayesPrism accepts linear TPM when raw counts are
# unavailable.

suppressPackageStartupMessages({
  library(data.table)
  library(immunedeconv)
})

# data.table evaluates column names within its own data mask.
utils::globalVariables(c("CELL_TYPE"))

# ------------------------------- Configuration -------------------------------

project_dir <- Sys.getenv(
  "CLL_TME_DIR",
  unset = path.expand("~/research/CLL_TME")
)
project_dir <- normalizePath(project_dir, mustWork = TRUE)
setwd(project_dir)

expression_file <- "CLLCMA_alpha/bulk/cll_broad_2022/data_mrna_seq_tpm.txt"
lm22_file <- "reference/LM22.txt"
cibersort_script <- "CIBERSORT.R"

# CLLCMA_alpha GSE-subsampled BayesPrism export (step_4).
sc_ref_dir <- "CLLCMA_alpha/integrated_annotated/bayesprism"
sc_counts_file <- file.path(sc_ref_dir, "sc_counts.rds")
sc_type_file <- file.path(sc_ref_dir, "cell_type_labels.rds")
sc_state_file <- file.path(sc_ref_dir, "cell_state_labels.rds")

output_dir <- "CLLCMA_alpha/results/deconvolution/cll_broad_2022"

# Remove methods from this vector when only a subset should be run.
methods_to_run <- c(
  "quantiseq",
  "xcell",
  "mcp_counter",
  "cibersort",
  "bayesprism"
)

n_cores <- 8L
overwrite_existing <- TRUE
save_full_bayesprism_object <- FALSE

# CIBERSORT: QN=FALSE for RNA-seq TPM; set perm >= 100 for empirical p-values.
cibersort_perm <- 0L
cibersort_qn <- FALSE
cibersort_markers_per_type <- 150L

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------ Helper functions -----------------------------

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "Package '", package, "' is required. Install it in the active R ",
      "environment and rerun the script.",
      call. = FALSE
    )
  }
}

as_sample_by_feature <- function(result) {
  result <- as.data.frame(result, check.names = FALSE)
  if (ncol(result) < 2L) {
    stop("Expected a cell-type column followed by sample columns.")
  }

  feature_names <- as.character(result[[1L]])
  values <- as.matrix(result[, -1L, drop = FALSE])
  storage.mode(values) <- "numeric"

  output <- t(values)
  colnames(output) <- feature_names
  rownames(output) <- colnames(result)[-1L]
  output
}

save_result <- function(result, method, estimate_type) {
  result <- as.matrix(result)
  storage.mode(result) <- "numeric"

  if (is.null(rownames(result)) || is.null(colnames(result))) {
    stop(method, " result must have sample row names and feature columns.")
  }
  if (anyDuplicated(rownames(result)) || anyDuplicated(colnames(result))) {
    stop(method, " result contains duplicated sample or feature names.")
  }

  wide <- as.data.table(result, keep.rownames = "SAMPLE_ID")
  wide_file <- file.path(
    output_dir,
    paste0(method, "_sample_by_feature.tsv")
  )
  fwrite(wide, wide_file, sep = "\t", na = "NA")
  saveRDS(
    result,
    file.path(output_dir, paste0(method, "_sample_by_feature.rds"))
  )

  long <- melt(
    wide,
    id.vars = "SAMPLE_ID",
    variable.name = "CELL_TYPE",
    value.name = "ESTIMATE"
  )
  long[, `:=`(
    METHOD = method,
    ESTIMATE_TYPE = estimate_type,
    CELL_TYPE = as.character(CELL_TYPE)
  )]
  setcolorder(
    long,
    c("SAMPLE_ID", "METHOD", "ESTIMATE_TYPE", "CELL_TYPE", "ESTIMATE")
  )
  fwrite(
    long,
    file.path(output_dir, paste0(method, "_long.tsv")),
    sep = "\t",
    na = "NA"
  )
  long
}

run_or_load <- function(method, estimate_type, runner) {
  tryCatch(
    {
      rds_file <- file.path(
        output_dir,
        paste0(method, "_sample_by_feature.rds")
      )

      if (file.exists(rds_file) && !overwrite_existing) {
        message("Loading existing ", method, " result: ", rds_file)
        result <- readRDS(rds_file)
        return(save_result(result, method, estimate_type))
      }

      message("Running ", method, " ...")
      result <- runner()
      save_result(result, method, estimate_type)
    },
    error = function(error) {
      error_message <- paste(method, "failed:", conditionMessage(error))
      warning(error_message, call. = FALSE)
      writeLines(
        error_message,
        file.path(output_dir, paste0(method, "_ERROR.txt"))
      )
      NULL
    }
  )
}

collapse_duplicate_genes <- function(input_file) {
  expression <- fread(input_file, na.strings = c("NA", ""))
  gene_column <- names(expression)[1L]

  expression <- expression[
    !is.na(get(gene_column)) & get(gene_column) != ""
  ]
  collapsed <- expression[
    ,
    lapply(.SD, mean, na.rm = TRUE),
    by = gene_column
  ]

  genes <- collapsed[[gene_column]]
  matrix <- as.matrix(collapsed[, -1L, with = FALSE])
  rownames(matrix) <- genes
  storage.mode(matrix) <- "numeric"

  if (any(!is.finite(matrix))) {
    stop("The collapsed expression matrix contains non-finite values.")
  }
  matrix
}

load_sc_reference <- function(collapse_rs_to_cll = FALSE) {
  require_package("Matrix")
  missing <- !file.exists(c(sc_counts_file, sc_type_file, sc_state_file))
  if (any(missing)) {
    stop(
      "Missing CLLCMA BayesPrism reference files under ", sc_ref_dir,
      ". Run CLLCMA_alpha/scripts/step_4_export_bayesprism_ref.R first."
    )
  }

  sc_counts <- readRDS(sc_counts_file)
  cell_type_labels <- as.character(readRDS(sc_type_file))
  cell_state_labels <- as.character(readRDS(sc_state_file))

  if (nrow(sc_counts) != length(cell_type_labels) ||
      length(cell_type_labels) != length(cell_state_labels)) {
    stop("sc_counts rows must match cell_type_labels and cell_state_labels.")
  }
  if (!"CLL" %in% cell_type_labels) {
    stop("Reference has no CLL cell_type; check step_4 export.")
  }

  # BayesPrism needs a single malignant key. Keep RS_* states but fold RS
  # into the CLL type so key = "CLL" covers all malignant reference cells.
  if (collapse_rs_to_cll && "RS" %in% cell_type_labels) {
    n_rs <- sum(cell_type_labels == "RS")
    cell_type_labels[cell_type_labels == "RS"] <- "CLL"
    message(
      "Collapsed ", n_rs, " RS reference cells into cell_type CLL ",
      "(cell_state labels unchanged)."
    )
  }

  list(
    counts = sc_counts,
    cell_type = cell_type_labels,
    cell_state = cell_state_labels
  )
}

# Build a CIBERSORT signature (genes x cell types) from scRNA mean CPM,
# selecting the top markers per type by fold-change vs other types.
build_sc_cibersort_signature <- function(
    sc_counts,
    cell_type_labels,
    n_markers_per_type = 150L
) {
  require_package("Matrix")

  types <- sort(unique(cell_type_labels))
  mean_counts <- vapply(
    types,
    function(tp) {
      idx <- which(cell_type_labels == tp)
      as.numeric(Matrix::colMeans(sc_counts[idx, , drop = FALSE]))
    },
    numeric(ncol(sc_counts))
  )
  rownames(mean_counts) <- colnames(sc_counts)
  colnames(mean_counts) <- types

  lib_sizes <- colSums(mean_counts)
  cpm <- sweep(mean_counts, 2L, lib_sizes, "/") * 1e6

  selected <- character()
  for (tp in types) {
    others <- setdiff(types, tp)
    other_mean <- rowMeans(cpm[, others, drop = FALSE])
    fold_change <- (cpm[, tp] + 1) / (other_mean + 1)
    top <- names(sort(fold_change, decreasing = TRUE))
    selected <- c(
      selected,
      head(top, n = min(n_markers_per_type, length(top)))
    )
  }
  selected <- unique(selected)
  signature <- cpm[selected, , drop = FALSE]
  storage.mode(signature) <- "numeric"
  signature
}

# ---------------------------- Expression matrix ------------------------------

# df is unlogged TPM with genes in rows and samples in columns. Raw counts are
# preferred for BayesPrism, but are not available in this cBioPortal download.
df <- collapse_duplicate_genes(expression_file)
message(
  "Expression matrix: ",
  nrow(df), " genes x ", ncol(df), " samples"
)

all_results <- list()

# -------------------------------- quanTIseq ----------------------------------

if ("quantiseq" %in% methods_to_run) {
  all_results[["quantiseq"]] <- run_or_load(
    method = "quantiseq",
    estimate_type = "fraction",
    runner = function() {
      as_sample_by_feature(
        immunedeconv::deconvolute(df, "quantiseq")
      )
    }
  )
}

# ---------------------------------- xCell ------------------------------------

if ("xcell" %in% methods_to_run) {
  require_package("xCell")
  require_package("GSVA")
  require_package("BiocParallel")

  # xCell 1.1.0 calls the GSVA API retired in GSVA >= 2.0. This compatibility
  # replacement keeps xCell's original algorithm while using ssgseaParam().
  if (utils::packageVersion("GSVA") >= "2.0.0") {
    raw_enrichment_analysis_fixed <- function(
        expr,
        signatures,
        genes,
        file.name = NULL,
        parallel.sz = 4,
        parallel.type = "SOCK"
    ) {
      shared_genes <- intersect(rownames(expr), genes)
      message("xCell shared genes: ", length(shared_genes))
      expr <- expr[shared_genes, , drop = FALSE]

      if (nrow(expr) < 5000L) {
        stop("Fewer than 5,000 genes overlap the xCell signatures.")
      }

      ranked_expr <- apply(expr, 2L, rank)
      params <- GSVA::ssgseaParam(
        exprData = ranked_expr,
        geneSets = signatures,
        normalize = FALSE
      )
      scores <- GSVA::gsva(
        params,
        verbose = FALSE,
        BPPARAM = BiocParallel::MulticoreParam(parallel.sz)
      )

      scores <- scores - apply(scores, 1L, min)
      cell_types <- unlist(strsplit(rownames(scores), "%", fixed = TRUE))
      cell_types <- cell_types[seq.int(1L, length(cell_types), by = 3L)]
      aggregated <- aggregate(scores ~ cell_types, FUN = mean)
      rownames(aggregated) <- aggregated[, 1L]
      aggregated[, -1L, drop = FALSE]
    }

    assignInNamespace(
      "rawEnrichmentAnalysis",
      raw_enrichment_analysis_fixed,
      ns = "xCell"
    )
  }

  all_results[["xcell"]] <- run_or_load(
    method = "xcell",
    estimate_type = "enrichment_score",
    runner = function() {
      as_sample_by_feature(
        immunedeconv::deconvolute(df, "xcell")
      )
    }
  )
}

# ------------------------------- MCP-counter ---------------------------------

if ("mcp_counter" %in% methods_to_run) {
  require_package("MCPcounter")

  all_results[["mcp_counter"]] <- run_or_load(
    method = "mcp_counter",
    estimate_type = "abundance_score",
    runner = function() {
      mcp_signature_file <- file.path(
        "reference",
        "MCPcounter_genes.txt"
      )
      if (!file.exists(mcp_signature_file)) {
        mcp_signature_url <- paste0(
          "https://raw.githubusercontent.com/",
          "ebecht/MCPcounter/master/Signatures/genes.txt"
        )
        message("Downloading MCP-counter signatures over HTTPS.")
        download.file(
          mcp_signature_url,
          destfile = mcp_signature_file,
          method = "libcurl",
          quiet = FALSE
        )
      }
      mcp_genes <- read.delim(
        mcp_signature_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )

      as_sample_by_feature(
        immunedeconv::deconvolute(
          df,
          "mcp_counter",
          genes = mcp_genes
        )
      )
    }
  )
}

# -------------------------------- CIBERSORT ----------------------------------

if ("cibersort" %in% methods_to_run) {
  require_package("e1071")
  require_package("preprocessCore")
  require_package("Matrix")

  all_results[["cibersort"]] <- run_or_load(
    method = "cibersort",
    estimate_type = "fraction",
    runner = function() {
      if (!file.exists(cibersort_script)) {
        stop("CIBERSORT script not found: ", cibersort_script)
      }

      source(cibersort_script, local = TRUE)
      # Keep CLL and RS as separate signature columns for CIBERSORT.
      sc_ref <- load_sc_reference(collapse_rs_to_cll = FALSE)

      signature <- build_sc_cibersort_signature(
        sc_counts = sc_ref$counts,
        cell_type_labels = sc_ref$cell_type,
        n_markers_per_type = cibersort_markers_per_type
      )

      common_genes <- intersect(rownames(df), rownames(signature))
      if (length(common_genes) < 50L) {
        stop(
          "Fewer than 50 signature genes overlap the expression matrix (",
          length(common_genes), ")."
        )
      }
      message(
        "CIBERSORT sc signature: ",
        ncol(signature), " cell types, ",
        length(common_genes), " / ", nrow(signature),
        " marker genes in bulk"
      )

      signature <- signature[common_genes, , drop = FALSE]
      sig_file <- file.path(output_dir, "cibersort_sc_signature.tsv")
      fwrite(
        as.data.table(signature, keep.rownames = "Gene symbol"),
        sig_file,
        sep = "\t"
      )

      mixture_file <- file.path(output_dir, "cibersort_mixture_tpm.tsv")
      mixture <- as.data.table(
        df[common_genes, , drop = FALSE],
        keep.rownames = "Gene symbol"
      )
      fwrite(mixture, mixture_file, sep = "\t")

      sig_abs <- normalizePath(sig_file, mustWork = TRUE)
      mixture_abs <- normalizePath(mixture_file, mustWork = TRUE)

      # CIBERSORT.R writes CIBERSORT-Results.txt to the working directory.
      old_wd <- getwd()
      on.exit(setwd(old_wd), add = TRUE)
      setwd(output_dir)

      cs <- CIBERSORT(
        sig_matrix = sig_abs,
        mixture_file = mixture_abs,
        perm = cibersort_perm,
        QN = cibersort_qn
      )

      setwd(old_wd)

      raw_results <- file.path(output_dir, "CIBERSORT-Results.txt")
      if (file.exists(raw_results)) {
        file.rename(
          raw_results,
          file.path(output_dir, "cibersort_raw_results.txt")
        )
      }

      qc_cols <- c("P-value", "Correlation", "RMSE")
      feature_cols <- setdiff(colnames(cs), qc_cols)
      fractions <- cs[, feature_cols, drop = FALSE]
      storage.mode(fractions) <- "numeric"

      qc <- cs[, intersect(qc_cols, colnames(cs)), drop = FALSE]
      if (ncol(qc) > 0L) {
        fwrite(
          as.data.table(qc, keep.rownames = "SAMPLE_ID"),
          file.path(output_dir, "cibersort_fit_metrics.tsv"),
          sep = "\t"
        )
      }

      fractions
    }
  )
}

# -------------------------------- BayesPrism ---------------------------------

if ("bayesprism" %in% methods_to_run) {
  require_package("BayesPrism")
  require_package("Matrix")

  all_results[["bayesprism"]] <- run_or_load(
    method = "bayesprism",
    estimate_type = "fraction",
    runner = function() {
      sc_ref <- load_sc_reference(collapse_rs_to_cll = TRUE)
      reference <- sc_ref$counts

      # BayesPrism: mixture is sample-by-gene; reference is cell-by-gene.
      mixture <- t(df)
      common_genes <- intersect(colnames(reference), colnames(mixture))
      if (length(common_genes) < 1000L) {
        stop(
          "Too few shared genes between sc reference and bulk (",
          length(common_genes), ")."
        )
      }
      message(
        "BayesPrism CLLCMA sc reference: ",
        nrow(reference), " cells, ",
        length(common_genes), " shared genes; key = CLL"
      )

      # Densify only the shared-gene subset to limit memory.
      reference <- as.matrix(reference[, common_genes, drop = FALSE])
      storage.mode(reference) <- "numeric"
      mixture <- mixture[, common_genes, drop = FALSE]

      cell_type_labels <- sc_ref$cell_type
      cell_state_labels <- sc_ref$cell_state

      stopifnot(
        identical(colnames(reference), colnames(mixture)),
        length(cell_type_labels) == nrow(reference),
        length(cell_state_labels) == nrow(reference),
        all(reference >= 0),
        all(mixture >= 0),
        "CLL" %in% cell_type_labels
      )

      prism <- BayesPrism::new.prism(
        reference = reference,
        mixture = mixture,
        input.type = "count.matrix",
        cell.type.labels = cell_type_labels,
        cell.state.labels = cell_state_labels,
        key = "CLL",
        outlier.cut = 0.01,
        outlier.fraction = 0.1
      )
      fit <- BayesPrism::run.prism(
        prism = prism,
        n.cores = n_cores
      )

      if (save_full_bayesprism_object) {
        saveRDS(
          fit,
          file.path(output_dir, "bayesprism_full_fit.rds"),
          compress = FALSE
        )
      }

      # Also export state-level fractions for patient-specific CLL clones.
      state_frac <- BayesPrism::get.fraction(
        bp = fit,
        which.theta = "final",
        state.or.type = "state"
      )
      fwrite(
        as.data.table(state_frac, keep.rownames = "SAMPLE_ID"),
        file.path(output_dir, "bayesprism_state_sample_by_feature.tsv"),
        sep = "\t"
      )
      saveRDS(
        state_frac,
        file.path(output_dir, "bayesprism_state_sample_by_feature.rds")
      )

      BayesPrism::get.fraction(
        bp = fit,
        which.theta = "final",
        state.or.type = "type"
      )
    }
  )
}

# ------------------------------- Combined data -------------------------------

if (!length(all_results)) {
  stop("No deconvolution methods were selected.")
}

combined_long <- rbindlist(all_results, use.names = TRUE)
setorder(combined_long, METHOD, CELL_TYPE, SAMPLE_ID)

fwrite(
  combined_long,
  file.path(output_dir, "all_methods_long.tsv"),
  sep = "\t",
  na = "NA"
)
saveRDS(
  combined_long,
  file.path(output_dir, "all_methods_long.rds")
)

method_summary <- combined_long[
  ,
  .(
    N_SAMPLES = uniqueN(SAMPLE_ID),
    N_CELL_TYPES = uniqueN(CELL_TYPE),
    N_MISSING = sum(is.na(ESTIMATE))
  ),
  by = .(METHOD, ESTIMATE_TYPE)
]
fwrite(
  method_summary,
  file.path(output_dir, "method_summary.tsv"),
  sep = "\t"
)

message("Finished. Results written to: ", normalizePath(output_dir))
