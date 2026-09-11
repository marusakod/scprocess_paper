suppressPackageStartupMessages({
  library("data.table")
  library("ggplot2")
  library("yaml")
  library("scales")
  library("magrittr")
  library("ggbeeswarm")
  library("stringr")
  library("ComplexHeatmap")
  library("circlize")
  library("SingleCellExperiment")
  library("edgeR")
  library("patchwork")
})

COMBOS = c("sample_excl_ambient", "sample_keep_ambient",
  "all_excl_ambient", "all_keep_ambient")

GSEA_REGEX = "^(HALLMARK_|GOBP_|GOCC_|GOMF_|BIOCARTA_|REACTOME_|KEGG_)(.+)"


# utils
cached_rds <- function(path, code, overwrite = FALSE, envir = parent.frame()) {
  if (!overwrite && file.exists(path)) {
    message("Loading cached: ", path)
    return(readRDS(path))
  }
  result = eval(substitute(code), envir = envir)
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(result, path)
  message("Saved cache: ", path)
  result
}

.stop_missing <- function(f, label) {
  if (!file.exists(f)) stop(sprintf("Required file not found (%s):\n  %s", label, f))
}

.parse_config <- function(config_f) {
  cfg = read_yaml(config_f)
  list(
    proj_dir   = cfg$project$proj_dir,
    short_tag  = cfg$project$short_tag,
    full_tag   = cfg$project$full_tag,
    date_stamp = cfg$project$date_stamp,
    tag        = paste0(cfg$project$full_tag, "_", cfg$project$date_stamp),
    proj_name  = basename(cfg$project$proj_dir)
  )
}

# data loading
load_hvg_comparison <- function(config_f, out_dir, zoom = NULL) {
  p = .parse_config(config_f)

  if (!is.null(zoom)) {
    prefix  = file.path(out_dir, p$proj_name, paste0("zoom_", zoom))
    edger_f = file.path(prefix, "edger_dt.csv.gz")
  } else {
    prefix  = file.path(out_dir, p$proj_name)
    edger_f = file.path(prefix, "edger_dt.csv.gz")
  }

  hvg_fs = setNames(file.path(prefix, COMBOS, "hvg_dt.csv.gz"), COMBOS)
  
  # check if all files exist
  for (nm in COMBOS) .stop_missing(hvg_fs[nm], paste("hvg_dt:", nm))
  .stop_missing(edger_f, "edger_empty_genes")
  
  # get highly variable genes for all COMBOS
  hvg_ls = lapply(COMBOS, function(nm) {
    fread(hvg_fs[nm])[highly_variable == TRUE, gene_id]
  }) %>% setNames(COMBOS)

  list(
    hvg_ls  = hvg_ls,
    edger_dt = fread(edger_f),
    meta     = list(full_tag = p$full_tag, proj_name = p$proj_name, zoom = zoom)
  )
}


load_pseudobulk <- function(gene_ids, cell_types, config_ls, out_dir, pseudo_count = 10, per_sample = FALSE){

  gene_ids     = unique(gene_ids)

  # download tarball with all pseudobulks from github repo
  url = "https://github.com/marusakod/scprocess_paper/releases/download/v0.1.0/pb_celltypes_data.tar.gz"
  temp_tar = tempfile(fileext = ".tar.gz")
  download.file(url, destfile = temp_tar, mode = "wb")

  message("Extracting files from pseudobulk tarball...")
  untar(temp_tar, exdir = ".")
  unlink(temp_tar)

  pb_dt = config_ls %>% lapply(function(config_f) {
    p = .parse_config(config_f)

    pb_f = file.path(out_dir, p$proj_name, "pb_celltypes.rds")
    if (!file.exists(pb_f)) {
      message("  [SKIP] missing pseudobulk celltypes file (", p$proj_name, "): ", pb_f)
      return(NULL)
    }
    pb = readRDS(pb_f)

    cell_types %>% lapply(function(ct){
      assay_nm = gsub(" ", "_", ct)
      if (!assay_nm %in% SummarizedExperiment::assayNames(pb)) {
        message("  [SKIP] missing cell type assay (", p$proj_name, ", ", ct, "): ", assay_nm)
        return(NULL)
      }

      mat = SummarizedExperiment::assay(pb, assay_nm)
      mat = mat[, colSums(mat) > 0, drop = FALSE]
      if (ncol(mat) == 0) return(NULL)

      present = intersect(gene_ids, rownames(mat))
      if (length(present) == 0) return(NULL)

      # library size = edgeR TMM-normalized effective library size
      dge      = suppressMessages(normLibSizes(DGEList(as.matrix(mat), remove.zeros = TRUE)))
      lib_size = getNormLibSizes(dge)

      sub_mat  = as.matrix(mat[present, , drop = FALSE])
      logcpm   = log(sweep(sub_mat, 2, lib_size, "/") * 1e6 + pseudo_count)

      dt = as.data.table(logcpm, keep.rownames = "gene_id") %>%
        melt(id.vars = "gene_id", variable.name = "sample_id", value.name = "logcpm")

      if (per_sample) {
        dt[, `:=`(cell_type = ct, study = p$proj_name)]
      } else {
        dt[, .(avg_logcpm = mean(logcpm)), by = gene_id] %>%
          .[, `:=`(cell_type = ct, study = p$proj_name)]
      }
     }) %>% rbindlist()
  }) %>% rbindlist()

  return(pb_dt)
}


load_fgsea_comparison <- function(config_f, out_dir, path_set = "go_bp") {
  p = .parse_config(config_f)

  fgsea_f = file.path(out_dir, p$proj_name, sprintf("fgsea_empties_%s.csv.gz", path_set))
  .stop_missing(fgsea_f, paste("fgsea_empties:", p$proj_name, path_set))

  dt = fread(fgsea_f) %>%
   .[, project := p$proj_name]

  return(dt)
}

# figures
.get_ambient_prop_dt <- function(dat_ls) {

  dt = lapply(dat_ls, function(dat) {
    ambient_gs = dat$edger_dt[is_ambient == TRUE, gene_id]
    tmp_dt = lapply(c('sample_keep_ambient', 'all_keep_ambient'), function(m){
    hvgs = dat$hvg_ls[[m]]
    ambient_hvgs = intersect(hvgs, ambient_gs)
    data.table(
      project   = dat$meta$proj_name,
      zoom      = dat$meta$zoom %||% "all cells",
      prop_ambient_hvgs = length(ambient_hvgs)/length(hvgs), 
      method    = m
    )
  }) %>% rbindlist()
  }) %>% rbindlist()
  
}

plot_prop_ambient <- function(dat_ls) {
  dt = .get_ambient_prop_dt(dat_ls) %>%
    .[, method := factor(method, levels = c("sample_keep_ambient", "all_keep_ambient"))] %>%
    .[zoom != "all cells"]

  stats_dt = copy(dt)[, .(
    med = median(prop_ambient_hvgs, na.rm = TRUE),
    q05 = quantile(prop_ambient_hvgs, 0.05, na.rm = TRUE),
    q95 = quantile(prop_ambient_hvgs, 0.95, na.rm = TRUE)
  ), by = .(method, zoom)]

  ggplot(dt, aes(x = method, y = prop_ambient_hvgs, fill = method)) +
    geom_quasirandom(shape = 21, size = 2, width = 0.2, alpha = 0.8) +
    theme_minimal() +
    facet_wrap(~zoom, nrow = 1, scales = "fixed", strip.position = "bottom") +
    scale_fill_brewer(palette = "Set2",
      labels = c(sample_keep_ambient = "batch-aware HVG selection", all_keep_ambient = "batch-unaware HVG selection")) +
    theme(
      axis.text.x      = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.placement  = "outside",
      strip.background = element_blank(),
      axis.line.y      = element_line(color = "black", linewidth = 0.5),
      axis.line.x      = element_line(color = "black", linewidth = 0.5),
      strip.text       = element_text(size = 11),
      axis.text.y      = element_text(size = 10),
      axis.title.y     = element_text(size = 11),
      legend.text      = element_text(size = 11),
      legend.position  = 'bottom'
    ) +
    scale_y_continuous(limits = c(0, 1), breaks = scales::pretty_breaks()) +
    labs(x = NULL, y = "Proportion of HVGs that are ambient", fill = NULL) +
    geom_segment(
      data = stats_dt,
      aes(x = method, xend = method, y = q05, yend = q95),
      linewidth = 0.5, colour = "black", inherit.aes = FALSE
    ) +
    geom_point(
      data = stats_dt,
      aes(x = method, y = med),
      colour = "black", size = 2, inherit.aes = FALSE
    )
}


.common_ambient_genes_dt <- function(dat_ls, n_top = 10, hvg_method = "all_keep_ambient") {
  flag_dt = dat_ls %>% lapply(function(dat) {
    # get ambient genes that are also hvgs
    hvg_ids = dat$hvg_ls[[hvg_method]]
    dat$edger_dt[is_ambient == TRUE & gene_id %chin% hvg_ids, .(gene_id, logFC)]
  }) %>% rbindlist()

  if (nrow(flag_dt) == 0) return(flag_dt[, .(gene_id, n_studies = integer(), mean_logFC = numeric())])
  
  # rank by how many studies flag them as ambient + hvg, with ties broken by average logFC across those studies
  flag_dt[, .(n_studies = .N, mean_logFC = mean(logFC)), by = gene_id] %>%
    setorder(-n_studies, -mean_logFC) %>%
    .[seq_len(min(n_top, .N))] %>%
    .[, gene_rank := seq_len(.N)]
}

# top common ambient + HVG genes per lineage
common_ambient_gene_table <- function(dat_ls, n_top = 20, hvg_method = c("sample_keep_ambient", "all_keep_ambient")){
  
  method = match.arg(hvg_method)
  zooms  = sapply(dat_ls, function(d) d$meta$zoom %||% "all cells")

  gene_dt = rbindlist(lapply(unique(zooms), function(z) {
    zoom_dat_ls = dat_ls[zooms == z]
    top_dt      = .common_ambient_genes_dt(zoom_dat_ls, n_top = n_top, hvg_method = hvg_method)
    if (nrow(top_dt) == 0) return(NULL)
    top_dt[, zoom := z]
  }))

  gene_dt[, symbol  := gene_id %>% str_extract(".+(?=_ENS)")]
  gene_dt[is.na(symbol), symbol := gene_id]

  return (gene_dt)
}


COMMON_AMBIENT_CELL_ORDER = c("myeloid", "oligos", "astrocytes", "neurons", "vascular", "all cells")
COMMON_AMBIENT_GLIAL_TYPES = c("myeloid", "oligos", "astrocytes")

.common_ambient_cell_levels <- function(vals) {
  c(intersect(COMMON_AMBIENT_CELL_ORDER, unique(vals)),
    setdiff(unique(vals), COMMON_AMBIENT_CELL_ORDER))
}

.common_ambient_plot_dt <- function(gene_dt, expr_dt, use_zscore) {
  studies    = unique(expr_dt$study)
  cell_types = unique(expr_dt$cell_type)

  grid_dt = CJ(gene_id = unique(gene_dt$gene_id), study = studies, cell_type = cell_types)
  grid_dt = merge(grid_dt, gene_dt[, .(gene_id, zoom, gene_rank, symbol)], by = "gene_id",
    allow.cartesian = TRUE)
  plot_dt = merge(grid_dt, expr_dt, by = c("gene_id", "cell_type", "study"), all.x = TRUE)

  gene_order = plot_dt[, .(zoom, gene_rank, symbol)] %>% unique %>%
    setorder(zoom, gene_rank) %>% .[, paste(zoom, symbol, sep = "|")]
  plot_dt[, gene_facet_id := factor(paste(zoom, symbol, sep = "|"), levels = rev(gene_order))]
  plot_dt[, zoom := factor(zoom, levels = .common_ambient_cell_levels(zoom))]
  plot_dt[, cell_type := factor(cell_type, levels = .common_ambient_cell_levels(cell_type))]

  if (use_zscore) {
    plot_dt[, value := as.numeric(scale(avg_logcpm)), by = .(gene_id, zoom)]
    fill_lab = "row z-score\n(log CPM)\n(counts after CellBender)"
  } else {
    plot_dt[, value := avg_logcpm]
    fill_lab = "mean log CPM"
  }

  list(plot_dt = plot_dt, fill_lab = fill_lab)
}


.term_membership_mat <- function(symbols, term_dt) {
  term_dt = as.data.table(term_dt)[, .(gene = list(gene)), by = term]
  mat = vapply(term_dt$gene, function(gs) {
    as.integer(symbols %in% gs)
  }, numeric(length(symbols)))
  matrix(as.integer(mat), nrow = length(symbols), dimnames = list(symbols, make.unique(term_dt$term)))
}


.build_common_ambient_complexheatmap <- function(plot_dt, use_zscore, fill_lab, term_dt = NULL) {
  row_ids = levels(plot_dt$gene_facet_id)

  col_dt = unique(plot_dt[, .(cell_type, study)]) %>% setorder(cell_type, study)
  col_dt[, col_id := paste(cell_type, study, sep = "|")]

  mat = matrix(NA_real_, nrow = length(row_ids), ncol = nrow(col_dt),
    dimnames = list(row_ids, col_dt$col_id))
  plot_dt[, col_id := paste(cell_type, study, sep = "|")]
  mat[cbind(as.character(plot_dt$gene_facet_id), plot_dt$col_id)] = plot_dt$value

  # the (gene, study, cell_type) grid built upstream is a full cross product, so it
  # includes study/cell_type combos where that zoom's pseudobulk was never computed at
  # all -- those columns/rows are NA for every entry. Drop them before display.

  keep_cols = colSums(!is.na(mat)) > 0
  keep_rows = rowSums(!is.na(mat)) > 0
  mat     = mat[keep_rows, keep_cols, drop = FALSE]
  col_dt  = col_dt[keep_cols]
  row_ids = row_ids[keep_rows]

  row_zoom  = sub("\\|.*", "", row_ids)
  row_split = factor(row_zoom, levels = .common_ambient_cell_levels(row_zoom))
  col_split = factor(col_dt$cell_type, levels = .common_ambient_cell_levels(col_dt$cell_type))

  pal = rev(RColorBrewer::brewer.pal(11, "RdBu"))
  if (use_zscore) {
    zlim = quantile(abs(mat), 0.99, na.rm = TRUE)
    col_fun = colorRamp2(seq(-zlim, zlim, length.out = length(pal)), pal)
  } else {
    rng = range(mat, na.rm = TRUE)
    col_fun = colorRamp2(seq(rng[1], rng[2], length.out = length(pal)), pal)
  }

  # hclust can't handle NA distances (missing zoom pseudobulk for a study/cell_type);
  # impute per-column mean for clustering purposes only, the displayed tile stays NA
  impute_for_clust = function(m) {
    m_imp = apply(m, 2, function(col) {
      col[is.na(col)] = mean(col, na.rm = TRUE)
      col
    })
    m_imp[is.nan(m_imp)] = 0
    m_imp
  }
  clust_rows_fun = function(m) stats::hclust(stats::dist(impute_for_clust(m)), method = "average")

  # outline the (zoom == cell_type) glial diagonal panels
  diagonal_layer_fun = function(j, i, x, y, w, h, fill, slice_r, slice_c) {
    row_grp = levels(row_split)[slice_r]
    col_grp = levels(col_split)[slice_c]
    if (row_grp %in% COMMON_AMBIENT_GLIAL_TYPES && identical(row_grp, col_grp)) {
      grid::grid.rect(gp = grid::gpar(fill = NA, col = "black", lwd = 3))
    }
  }

  right_annot = NULL
  top_pad = unit(2, "mm")
  if (!is.null(term_dt) && nrow(term_dt) > 0) {
    row_symbols = sub(".*\\|", "", row_ids)
    term_mat = .term_membership_mat(row_symbols, term_dt)
    term_df = as.data.frame(lapply(as.data.frame(term_mat), factor, levels = c(0, 1)))
    colnames(term_df) = colnames(term_mat)
    annot_fontsize = 9
    right_annot = rowAnnotation(
      df = term_df,
      col = setNames(rep(list(c(`0` = "white", `1` = "black")), ncol(term_df)), colnames(term_df)),
      show_legend         = FALSE,
      annotation_name_side = "top",
      annotation_name_rot = 90,
      annotation_name_gp  = gpar(fontsize = annot_fontsize),
      simple_anno_size    = unit(4, "mm"),
      gp                  = gpar(col = "black", lwd = 0.5)
    )
    name_widths_mm = vapply(colnames(term_df), function(s) {
      grid::convertWidth(grid::grobWidth(grid::textGrob(s, gp = gpar(fontsize = annot_fontsize))),
        "mm", valueOnly = TRUE)
    }, numeric(1))
    top_pad = unit(max(name_widths_mm) + 6, "mm")
  }

  ht = Heatmap(mat,
    name              = fill_lab,
    col               = col_fun,
    na_col            = "grey85",
    row_split         = row_split,
    cluster_row_slices = FALSE,
    cluster_rows      = clust_rows_fun,
    show_row_dend     = TRUE,
    row_labels        = sub(".*\\|", "", row_ids),
    row_names_gp      = gpar(fontsize = 8),
    row_title_rot     = 0,
    row_title_gp      = gpar(fontsize = 10, fontface = "bold"),
    column_split      = col_split,
    cluster_columns   = FALSE,
    column_labels     = col_dt$study,
    column_names_gp   = gpar(fontsize = 8),
    column_title_gp   = gpar(fontsize = 10, fontface = "bold"),
    border            = TRUE,
    rect_gp           = gpar(col = "white", lwd = 0.5),
    layer_fun         = diagonal_layer_fun,
    right_annotation  = right_annot,
    heatmap_legend_param = list(title = fill_lab)
  )

  ComplexHeatmap::draw(ht,
    column_title      = "Columns are cell types for gene expression",
    column_title_side = "top",
    row_title         = "Rows are cell types for gene selection",
    row_title_side    = "left",
    padding           = unit.c(unit(2, "mm"), unit(2, "mm"), top_pad, unit(2, "mm"))
  )
  invisible(NULL)
}


plot_common_ambient_heatmap <- function(gene_dt, expr_dt, use_zscore = FALSE, term_dt = NULL) {

  prep = .common_ambient_plot_dt(gene_dt, expr_dt, use_zscore)
  .build_common_ambient_complexheatmap(prep$plot_dt, use_zscore, prep$fill_lab, term_dt = term_dt)
}


plot_gsea_dotplot_by_lineage <- function(dat_ls, sel_lineage, n_top_up = 5,
  n_top_down = 5, gsea_cut = 0.05, maxp_cut = 0.1, n_chars = 50,
  max_nes = Inf, min_log10_padj = -10, project_order = NULL) {

  gsea_dt = rbindlist(dat_ls, fill = TRUE) %>%
    .[lineage == sel_lineage]
  
  if (nrow(gsea_dt) == 0) return(NULL)

  base_dt = copy(gsea_dt)[main_path == TRUE] %>%
    .[, min_p := min(padj, na.rm = TRUE), by = pathway] %>%
    .[min_p < maxp_cut]
  if (nrow(base_dt) == 0) return(NULL)

  path_summary = base_dt[, .(min_padj = min(padj, na.rm = TRUE),
    med_nes = median(NES, na.rm = TRUE)), by = pathway]
  top_up = path_summary[med_nes > 0][order(min_padj)][seq_len(min(.N, n_top_up))]
  top_dn = path_summary[med_nes < 0][order(min_padj)][seq_len(min(.N, n_top_down))]
  top_paths = unique(c(top_up$pathway, top_dn$pathway))

  .clean_pathway = function(p) str_match(p, GSEA_REGEX)[, 3] %>%
    tolower %>% str_replace_all("_", " ") %>% str_sub(1, n_chars)

  plot_dt = gsea_dt[pathway %in% top_paths] %>%
    .[is.na(padj), padj := 1] %>%
    .[, path_short := .clean_pathway(pathway)] %>%
    .[, pathway    := str_match(pathway, GSEA_REGEX)[, 3] %>%
      tolower %>% str_replace_all("_", " ")] %>%
    .[, signif     := ifelse(padj < gsea_cut, "significant", "not")]

  if (!is.null(project_order))
    plot_dt[, project := factor(project, levels = project_order)]

  order_dt = plot_dt[, .(med_nes = median(NES, na.rm = TRUE)), by = path_short] %>%
    setorder(med_nes)
  plot_dt[, path_short := factor(path_short, levels = order_dt$path_short)]

  res     = 0.5
  max_nes = min(ceiling(max(abs(plot_dt$NES)) * res) / res, max_nes)
  plot_dt[, nes_trunc        := sign(NES) * pmin(abs(NES), max_nes)]
  plot_dt[, log10_padj_trunc := pmax(log10(padj), min_log10_padj)]
  plot_dt[order(NES)]

  ggplot(plot_dt) +
    aes(x = project, y = path_short, fill = nes_trunc,
        size = -log10_padj_trunc, alpha = signif) +
    geom_point(shape = 21, colour = "black") +
    scale_fill_distiller(palette = "RdBu", limits = c(-max_nes, max_nes),
      breaks = pretty_breaks()) +
    scale_alpha_manual(values = c(significant = 1, not = 0.5)) +
    scale_size(range = c(1, 6), breaks = pretty_breaks()) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 13),
          panel.grid  = element_blank()) +
    labs(x = "Project", y = NULL, title = sel_lineage,
         fill  = "NES",
         size  = "-log10(padj)",
         alpha = sprintf("Signif.\n(%d%%)", round(100 * gsea_cut)))
}


leading_edge_overlap_dt <- function(fgsea_data, target_genes, gsea_cut = 0.05) {
  gsea_dt = lapply(fgsea_data, FUN = rbindlist, fill = TRUE) %>% 
    rbindlist(fill = TRUE)

  sig_dt = gsea_dt[main_path == TRUE & padj < gsea_cut] %>%
  .[, leading_edge := strsplit(leadingEdge, "\\|")]
    
  target_genes = unique(target_genes)
  pooled_dt = sig_dt[, .(
    n_rows       = .N,                               
    leading_edge = list(unique(unlist(leading_edge)))
  ), by = pathway] %>%
  .[, n_leading_edge := lengths(leading_edge)] %>%
  .[, matched       := lapply(leading_edge, function(g) intersect(g, target_genes))] %>%
  .[, n_matched     := lengths(matched)] %>%
  .[, matched_genes := sapply(matched, paste, collapse = "/")] %>%
  .[, coverage      := round((n_matched / length(target_genes)), 2)] %>%  # fraction of target_genes hit
  .[, precision     := round((n_matched / n_leading_edge), 3)] %>%        # fraction of leading edge that's target_genes
  .[, term := str_match(pathway, GSEA_REGEX)[, 3] %>% tolower %>% str_replace_all("_", " ")] %>%
  .[n_matched > 0] %>%
  .[order(-n_matched, -coverage)]
  
  return(pooled_dt)
}


greedy_nonredundant_terms <- function(leading_edge_dt, min_new = 2) {
  dt = copy(leading_edge_dt)
  setorder(dt, -n_matched, -coverage)
  dt[, matched_set := strsplit(matched_genes, "/")]

  covered = character(0)
  keep_rows = list()
  for (i in seq_len(nrow(dt))) {
    g = dt$matched_set[[i]]
    new_g = setdiff(g, covered)
    if (length(new_g) < min_new) next
    covered = union(covered, g)
    keep_rows[[length(keep_rows) + 1]] = data.table(row = i, n_new = length(new_g))
  }

  kept = rbindlist(keep_rows)
  out = dt[kept$row] %>%
  .[, n_new := kept$n_new] %>%
  .[, matched_set := NULL]
  
  return(out)
}