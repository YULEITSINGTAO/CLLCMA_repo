#!/usr/bin/env Rscript
# Manhattan + QQ plots for GEMMA or PLINK CIBERSORT GWAS results using topr.
#
# Usage (from phs000879_GWAS/):
#   Rscript plot_single_manhattan.R CD4_T
#   Rscript plot_single_manhattan.R all gemma_bayesprism_covar
#   Rscript plot_single_manhattan.R CLL gemma_bayesprism_covar 1e-5
#   Rscript plot_single_manhattan.R CD4_T plink_cibersort_linear 1e-5
#
# Args:
#   1) phenotype name or "all"
#   2) result directory (default: gemma_bayesprism_covar)
#   3) annotate p-value threshold (default: 1e-5)
#
# Also writes *_manhattan_compressed.png: piecewise y-axis that squeezes
# -log10(P) in [0, compress_cut] and expands the region above it.
#
# Supports:
#   - GEMMA:  output/assoc_<TRAIT>.assoc.txt   (chr, ps, p_wald)
#   - PLINK:  output/assoc_<TRAIT>.assoc.linear (CHR, BP, P; TEST==ADD)

suppressPackageStartupMessages({
  library(data.table)
  library(topr)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
pheno <- if (length(args) >= 1) args[[1]] else "CD4_T"
outdir <- if (length(args) >= 2) args[[2]] else "gemma_bayesprism_covar"
annotate_p <- if (length(args) >= 3) as.numeric(args[[3]]) else 1e-5
build <- 38L
# Compress background cloud below this -log10(P); expand hits above it.
compress_cut <- 5
compress_ratio <- 0.15  # [0,5] maps to 25% of an equal-scale interval
ymax_plot <- 8

assoc_dir <- file.path(outdir, "output")
man_dir <- file.path(outdir, "manhattan_plots")
qq_dir <- file.path(outdir, "qq_plots")
dir.create(man_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(qq_dir, showWarnings = FALSE, recursive = TRUE)

method_label <- if (grepl("plink", outdir, ignore.case = TRUE)) {
  "PLINK linear"
} else {
  "GEMMA LMM"
}
pheno_label <- if (grepl("bayesprism", outdir, ignore.case = TRUE)) {
  "BayesPrism INT · sex+age10+IGHV"
} else if (grepl("cibersort", outdir, ignore.case = TRUE)) {
  "CIBERSORT INT"
} else {
  outdir
}

default_phenos <- c("CLL", "CD4_T", "CD8_T", "DC", "Monocyte", "NK", "Progenitor", "RS")
pheno_list_file <- file.path(outdir, "phenotypes_to_run.txt")
if (file.exists(pheno_list_file)) {
  default_phenos <- scan(pheno_list_file, what = character(), quiet = TRUE)
}
phenos <- if (identical(pheno, "all")) default_phenos else pheno

cat("outdir=", outdir, "\n", sep = "")
cat("method=", method_label, "\n", sep = "")
cat("phenotypes=", paste(phenos, collapse = ", "), "\n", sep = "")
cat("annotate_p=", annotate_p, " build=GRCh", build, "\n", sep = "")

# For very large assoc files, keep all SNPs below a loose plot threshold
# plus a random subset of the rest so topr stays responsive.
downsample_for_plot <- function(dt, keep_p = 1e-3, n_bg = 800000L) {
  hit <- dt[P <= keep_p]
  bg <- dt[P > keep_p]
  if (nrow(bg) > n_bg) {
    set.seed(1L)
    bg <- bg[sample.int(nrow(bg), n_bg)]
  }
  rbindlist(list(hit, bg), use.names = TRUE)
}

# Piecewise y transform: squeeze [0, cut], keep linear spacing above cut.
# Example with cut=5, ratio=0.25:
#   data 0→0, 5→1.25, 10→6.25  (hits get ~4x more vertical room per unit)
compress_logp_trans <- function(cut = 5, ratio = 0.25) {
  force(cut)
  force(ratio)
  scales::trans_new(
    name = paste0("compress_logp_", cut),
    transform = function(y) {
      y <- as.numeric(y)
      ifelse(
        !is.finite(y), y,
        ifelse(y <= cut, y * ratio, cut * ratio + (y - cut))
      )
    },
    inverse = function(z) {
      z <- as.numeric(z)
      thresh <- cut * ratio
      ifelse(
        !is.finite(z), z,
        ifelse(z <= thresh, z / ratio, cut + (z - thresh))
      )
    },
    domain = c(0, Inf)
  )
}

apply_compressed_yaxis <- function(p, cut = compress_cut, ratio = compress_ratio,
                                   ymax = ymax_plot) {
  breaks <- unique(c(0, 2.5, cut, seq(cut + 1, ymax, by = 1)))
  breaks <- breaks[breaks <= ymax]
  p +
    ggplot2::scale_y_continuous(
      trans = compress_logp_trans(cut = cut, ratio = ratio),
      breaks = breaks,
      labels = as.character(breaks),
      limits = c(0, ymax),
      expand = ggplot2::expansion(mult = c(0.01, 0.05)),
      name = expression(-log[10](P))
    ) +
    ggplot2::geom_hline(
      yintercept = cut,
      linetype = "dotted",
      color = "grey30",
      linewidth = 0.4
    ) +
    ggplot2::labs(
      subtitle = paste0(
        "Compressed axis: -log10(P) < ", cut,
        " squeezed (ratio=", ratio, "); focus on SNPs ≥ ", cut
      )
    )
}

find_assoc_file <- function(feature) {
  candidates <- c(
    file.path(assoc_dir, paste0("assoc_", feature, ".assoc.linear")),
    file.path(assoc_dir, paste0("assoc_", feature, ".assoc.txt"))
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

read_assoc_for_plot <- function(f) {
  results <- fread(f)
  cols <- names(results)
  # PLINK linear
  if (all(c("CHR", "BP", "P") %in% cols)) {
    if ("TEST" %in% cols) {
      results <- results[TEST == "ADD"]
    }
    plot_data <- data.table(
      CHROM = as.integer(results$CHR),
      POS = as.numeric(results$BP),
      P = as.numeric(results$P)
    )
  } else if (all(c("chr", "ps", "p_wald") %in% cols)) {
    # GEMMA
    plot_data <- data.table(
      CHROM = as.integer(results$chr),
      POS = as.numeric(results$ps),
      P = as.numeric(results$p_wald)
    )
  } else {
    stop("Unrecognized association columns in ", f, ": ", paste(cols, collapse = ", "))
  }
  plot_data[
    !is.na(P) & P > 0 & !is.na(CHROM) & !is.na(POS) &
      CHROM >= 1L & CHROM <= 22L
  ]
}

plot_one <- function(feature) {
  f <- find_assoc_file(feature)
  if (is.na(f)) {
    warning("Missing association file for: ", feature)
    return(invisible(NULL))
  }
  cat("Plotting ", feature, " <- ", f, "\n", sep = "")

  plot_data <- read_assoc_for_plot(f)

  n_sig <- plot_data[P <= annotate_p, .N]
  cat("  SNPs=", nrow(plot_data), "  P<=", annotate_p, ": ", n_sig, "\n", sep = "")

  plot_df <- as.data.frame(downsample_for_plot(plot_data))
  cat("  plotted points after downsample=", nrow(plot_df), "\n", sep = "")

  colors <- get_topr_colors()
  title_extra <- paste0(method_label, " · ", pheno_label)

  man_png <- file.path(man_dir, paste0(feature, "_manhattan_topr.png"))
  png(man_png, width = 1800, height = 600, res = 150)
  p <- tryCatch(
    {
      manhattanExtra(
        plot_df,
        genome_wide_thresh = 5e-8,
        suggestive_thresh = annotate_p,
        flank_size = 1e6,
        region_size = 1e6,
        annotate = annotate_p,
        color = colors,
        build = build,
        title = paste0("Manhattan - ", feature, " (", title_extra, ")"),
        protein_coding_only = TRUE,
        annotate_with = "Gene_Symbol",
        max.overlaps = 30,
        label_size = 3.2,
        size = 1.4,
        alpha = 0.85,
        ymax = ymax_plot
      )
    },
    error = function(e) {
      message("manhattanExtra failed (", conditionMessage(e), "); falling back to manhattan()")
      manhattan(
        plot_df,
        annotate = annotate_p,
        color = colors,
        build = build,
        title = paste0("Manhattan - ", feature, " (", title_extra, ")"),
        sign_thresh = c(5e-8, annotate_p),
        sign_thresh_color = c("red", "blue"),
        protein_coding_only = TRUE,
        annotate_with = "Gene_Symbol",
        region_size = 1e6,
        max.overlaps = 30,
        label_size = 3.2,
        size = 1.4,
        alpha = 0.85,
        ymax = ymax_plot
      )
    }
  )
  print(p)
  dev.off()
  cat("  wrote ", man_png, "\n", sep = "")

  man_comp_png <- file.path(man_dir, paste0(feature, "_manhattan_compressed.png"))
  png(man_comp_png, width = 1800, height = 600, res = 150)
  print(apply_compressed_yaxis(p))
  dev.off()
  cat("  wrote ", man_comp_png, "\n", sep = "")

  man2_png <- file.path(man_dir, paste0(feature, "_manhattan.png"))
  png(man2_png, width = 1800, height = 600, res = 150)
  p2 <- manhattan(
    plot_df,
    annotate = annotate_p,
    color = colors,
    build = build,
    title = paste0(feature, " ", title_extra, " (annotate P<", annotate_p, ")"),
    sign_thresh = c(5e-8, annotate_p),
    sign_thresh_color = c("firebrick", "steelblue"),
    protein_coding_only = TRUE,
    annotate_with = "Gene_Symbol",
    region_size = 1e6,
    max.overlaps = 40,
    nudge_y = 0.35,
    label_size = 3.0,
    size = 1.3,
    alpha = 0.9,
    ymax = ymax_plot
  )
  print(p2)
  dev.off()
  cat("  wrote ", man2_png, "\n", sep = "")

  man2_comp_png <- file.path(man_dir, paste0(feature, "_manhattan_compressed_basic.png"))
  png(man2_comp_png, width = 1800, height = 600, res = 150)
  print(apply_compressed_yaxis(p2))
  dev.off()
  cat("  wrote ", man2_comp_png, "\n", sep = "")

  qq_png <- file.path(qq_dir, paste0(feature, "_qq.png"))
  png(qq_png, width = 900, height = 600, res = 150)
  if ("qqplot" %in% getNamespaceExports("topr")) {
    print(qqplot(plot_data$P, title = paste0("QQ - ", feature, " (", method_label, ")")))
  } else {
    obs <- -log10(sort(plot_data$P))
    exp <- -log10(ppoints(length(obs)))
    plot(
      exp, obs,
      pch = 19, cex = 0.4, col = "grey25",
      xlab = expression(Expected ~ -log[10](P)),
      ylab = expression(Observed ~ -log[10](P)),
      main = paste0("QQ - ", feature, " (", method_label, ")")
    )
    abline(0, 1, col = "firebrick", lwd = 2)
  }
  dev.off()
  cat("  wrote ", qq_png, "\n", sep = "")

  leads <- plot_data[P <= annotate_p][order(P)]
  if (nrow(leads) > 0) {
    leads[, feature := feature]
    lead_file <- file.path(man_dir, paste0(feature, "_lead_snps.tsv"))
    fwrite(leads, lead_file, sep = "\t")
    cat("  wrote ", lead_file, " (", nrow(leads), " SNPs)\n", sep = "")
  }

  invisible(NULL)
}

for (p in phenos) {
  plot_one(p)
}
cat("Completed!\n")
