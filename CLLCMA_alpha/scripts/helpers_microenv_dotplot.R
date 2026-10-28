# Shared dot plot: mutation/alteration rate vs mean cell-type fraction.
#
# y = 100 * N_altered / (N_altered + N_WT)   (percent)
# x = 100 * mean CIBERSORT fraction among altered samples (percent)
# Per cell type: linear regression of mutation rate (%) ~ mean fraction (%)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

pct_axis_lab <- function(x) {
  paste0(format(x, trim = TRUE, scientific = FALSE, digits = 3), "%")
}

fmt_p <- function(p) {
  ifelse(
    is.na(p),
    "P = NA",
    ifelse(p < 0.001, sprintf("P = %.1e", p), sprintf("P = %.3f", p))
  )
}

fit_rate_vs_mean <- function(plot_dt) {
  plot_dt[
    ,
    {
      ok <- is.finite(MEAN_CT_PCT) & is.finite(MUT_RATE_PCT)
      x <- MEAN_CT_PCT[ok]
      y <- MUT_RATE_PCT[ok]
      n <- length(x)
      if (n < 5L || sd(x) == 0 || sd(y) == 0) {
        .(
          N = n,
          SLOPE = NA_real_,
          INTERCEPT = NA_real_,
          R2 = NA_real_,
          PEARSON_R = NA_real_,
          P_VALUE = NA_real_
        )
      } else {
        fit <- stats::lm(y ~ x)
        s <- summary(fit)
        ct <- stats::cor.test(x, y, method = "pearson")
        .(
          N = n,
          SLOPE = unname(stats::coef(fit)[[2L]]),
          INTERCEPT = unname(stats::coef(fit)[[1L]]),
          R2 = unname(s$r.squared),
          PEARSON_R = unname(ct$estimate),
          P_VALUE = unname(ct$p.value)
        )
      }
    },
    by = CELL_TYPE
  ]
}

plot_rate_vs_mean_celltype <- function(
  dt,
  gene_col,
  n_alt_col,
  n_wt_col,
  mean_alt_col,
  p_col,
  fdr_col,
  out_file,
  title,
  subtitle,
  highlight_genes = character(),
  cell_levels = character()
) {
  plot_dt <- data.table::copy(dt)
  plot_dt[, GENE := as.character(get(gene_col))]
  plot_dt[, MUT_RATE := get(n_alt_col) / (get(n_alt_col) + get(n_wt_col))]
  plot_dt[, MEAN_CT := get(mean_alt_col)]
  plot_dt[, PVAL := get(p_col)]
  plot_dt[, FDRVAL := get(fdr_col)]
  plot_dt[, N_ALT := get(n_alt_col)]
  plot_dt <- plot_dt[is.finite(MUT_RATE) & is.finite(MEAN_CT)]
  if (!nrow(plot_dt)) {
    warning("No rows for ", out_file)
    return(invisible(NULL))
  }

  plot_dt[, MUT_RATE_PCT := 100 * MUT_RATE]
  plot_dt[, MEAN_CT_PCT := 100 * MEAN_CT]

  cts <- unique(as.character(plot_dt$CELL_TYPE))
  lv <- c(intersect(cell_levels, cts), setdiff(cts, cell_levels))
  plot_dt[, CELL_TYPE := factor(CELL_TYPE, levels = lv)]
  plot_dt[, SIG := !is.na(FDRVAL) & FDRVAL < 0.05]
  # Label FDR hits plus a small driver set; cap at 6 per facet.
  plot_dt[
    ,
    SHOW_LABEL := {
      sig_genes <- unique(GENE[SIG])
      hi <- unique(GENE[GENE %in% highlight_genes])
      hi <- hi[!hi %in% sig_genes]
      if (length(sig_genes) + length(hi) > 6L) {
        n_hi <- max(0L, 6L - length(sig_genes))
        hi_ord <- order(-abs(CLIFFS_DELTA[match(hi, GENE)]))
        hi <- hi[hi_ord][seq_len(min(n_hi, length(hi)))]
      }
      GENE %in% c(sig_genes, hi)
    },
    by = CELL_TYPE
  ]
  plot_dt[, LAB := fifelse(SHOW_LABEL, GENE, NA_character_)]

  reg <- fit_rate_vs_mean(plot_dt)
  reg[, FDR := p.adjust(P_VALUE, method = "BH")]
  setorder(reg, P_VALUE)
  reg_file <- sub("\\.pdf$", "_regression.tsv", out_file)
  fwrite(reg, reg_file, sep = "\t")

  ann <- merge(
    reg,
    plot_dt[, .(
      X = min(MEAN_CT_PCT, na.rm = TRUE),
      Y = max(MUT_RATE_PCT, na.rm = TRUE)
    ), by = CELL_TYPE],
    by = "CELL_TYPE"
  )
  ann[, REG_LAB := sprintf(
    "R = %.2f\n%s",
    PEARSON_R,
    fmt_p(P_VALUE)
  )]
  ann[is.na(PEARSON_R), REG_LAB := ""]

  p <- ggplot(
    plot_dt,
    aes(x = MEAN_CT_PCT, y = MUT_RATE_PCT)
  ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      color = "grey25",
      fill = "grey75",
      linewidth = 0.45,
      alpha = 0.25
    ) +
    geom_point(
      aes(color = CLIFFS_DELTA, size = N_ALT, shape = SIG),
      alpha = 0.85
    ) +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17), guide = "none") +
    ggrepel::geom_text_repel(
      aes(label = LAB),
      size = 2.2,
      max.overlaps = 12,
      min.segment.length = 0,
      box.padding = 0.15,
      show.legend = FALSE,
      color = "grey20",
      na.rm = TRUE
    ) +
    geom_text(
      data = ann,
      aes(x = X, y = Y, label = REG_LAB),
      hjust = 0,
      vjust = 1,
      size = 2.3,
      color = "grey20",
      lineheight = 0.95,
      inherit.aes = FALSE
    ) +
    facet_wrap(~CELL_TYPE, scales = "free_x", ncol = 4) +
    scale_color_gradient2(
      low = "#2166AC", mid = "#BDBDBD", high = "#B2182B",
      midpoint = 0, limits = c(-0.6, 0.6),
      name = "Cliff's delta\n(mut vs WT)"
    ) +
    scale_size_continuous(range = c(1.4, 5.5), name = "N mutated") +
    scale_x_continuous(labels = pct_axis_lab) +
    scale_y_continuous(labels = pct_axis_lab, expand = expansion(mult = c(0.05, 0.18))) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Mean cell-type percentage (mutated samples)",
      y = "Mutation rate (%)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold", size = 8),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
  n_facet <- uniqueN(plot_dt$CELL_TYPE)
  ggsave(
    out_file,
    p,
    width = 13,
    height = max(7, 2.15 * ceiling(n_facet / 4) + 1.8)
  )
  invisible(list(plot = p, regression = reg))
}
