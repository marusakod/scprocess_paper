# Build per-project pseudobulk SingleCellExperiment objects from scprocess
# integration output, with one assay per cell type label

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(zellkonverter)
  library(data.table)
  library(BiocParallel)
  library(yaml)
  library(magrittr)
  library(stringr)
  library(Matrix)
  library(forcats)
  library(assertthat)
  library(purrr)
})

setDTthreads(1)

.parse_config <- function(config_f) {
  cfg = read_yaml(config_f)
  list(
    proj_dir   = cfg$project$proj_dir,
    short_tag  = cfg$project$short_tag,
    full_tag   = cfg$project$full_tag,
    date_stamp = cfg$project$date_stamp,
    tag        = paste0(cfg$project$full_tag, "_", cfg$project$date_stamp),
    proj_name  = basename(cfg$project$proj_dir),
    zoom       = cfg$zoom
  )
}

.get_zoom_names_and_lbls <- function(zoom_cfg_list) {
  lbls_dt = zoom_cfg_list %>% lapply(function(cfg_f){
    cfg = read_yaml(cfg_f)
    data.table(label = cfg$zoom$sel_labels, zoom = cfg$zoom$name)
  }) %>% rbindlist()
  
  return(lbls_dt)
}


aggregateData_datatable <- function(sce, by_vars = c("cluster", "sample_id"),
  fun = c("sum", "mean", "median", "prop.detected", "num.detected"), all_cls) {
  fun       = match.arg(fun)
  assy      = "X"

  counts    = assay(sce, "X") %>% as("TsparseMatrix")
  mat_dt    = data.table(
    i         = counts@i + 1,
    j         = counts@j + 1,
    count     = counts@x
    ) %>%
    .[, cell_id := colnames(sce)[j] ] %>%
    setkey("cell_id")

  by_1      = by_vars[[1]]
  by_2      = by_vars[[2]]
  by_dt     = data.table(
    cell_id   = colnames(sce),
    cluster   = factor(sce[[ by_1 ]]),
    sample_id = factor(sce[[ by_2 ]]) %>% fct_drop
    ) %>% setkey("cell_id")

  pb_dt     = mat_dt %>% merge(by_dt, by = "cell_id") %>%
    .[, .(sum = sum(count), .N), by = c("i", by_vars) ]

  n_cells_dt  = by_dt[, .(n_cells = .N), by = by_vars ]
  rows_dt     = data.table(i = seq.int(nrow(sce)), gene_id = rownames(sce))
  pb_dt       = pb_dt %>%
    merge(n_cells_dt, by = by_vars) %>%
    merge(rows_dt, by = "i") %>%
    .[, i             := NULL ] %>%
    .[, prop.detected := N / n_cells ]

  genes_ls    = rownames(sce)
  cl_ls       = levels(by_dt$cluster)
  samples     = levels(by_dt$sample_id)
  n_genes     = length(genes_ls)
  n_samples   = length(samples)

  mat_ls      = lapply(cl_ls, function(cl) {
    outs_mat    = pb_dt[ (cluster == cl) & (sample_id != "") ] %>%
      dcast.data.table( gene_id ~ sample_id, value.var = fun, fill = 0 ) %>%
      as.matrix(rownames = "gene_id")

    mat         = matrix(0, nrow=n_genes, ncol=n_samples) %>%
        set_rownames(genes_ls) %>%
        set_colnames(samples)
    mat[ rownames(outs_mat), colnames(outs_mat) ] = outs_mat

    return(mat)
    }) %>% setNames(cl_ls)

  md          = metadata(sce)
  md$agg_pars = list(assay = assy, by = by_vars, fun = fun)
  pb          = SingleCellExperiment(mat_ls,
    rowData = rowData(sce), metadata = md)
  cd          = data.frame(colData(sce)[, by_vars])

  cd[[by_1]]  = factor(cd[[by_1]], levels = all_cls)
  ns          = table(cd)
  if (length(by_vars) == 2) {
    ns    = asplit(ns, 2)
    ns    = purrr::map(ns, ~c(unclass(.)))
  } else {
    ns     = c(unclass(ns))
  }
  int_colData(pb)$n_cells = ns

  assert_that( all(rownames(pb) == rownames(sce)) )
  return(pb)
}


.make_one_pseudobulk_celltypes <- function(sel_sample, h5ad_paths, labels_dt, labels_with_zoom,
  agg_fn = aggregateData_datatable){
  
  message(sel_sample)
  h5ad_f  = h5ad_paths[[sel_sample]]
  #check if file is empty for that sample (sample has been excluded)
  if(file.size(h5ad_f) == 0){
     message("  sample ", sel_sample, " has been excluded -- skipping")
    return(NULL)
  }
  tmp_sce = readH5AD(h5ad_f)
  
  zoom_lbls      = zoom_lbls = unique(labels_dt$zoom)
  smpl_labels_dt = labels_dt[sample_id == sel_sample] %>% setkey(cell_id)
  keep_ids       = intersect(colnames(tmp_sce), smpl_labels_dt$cell_id)
  if (length(keep_ids) == 0) {
    message("  no confidently-labelled cells for ", sel_sample, " -- skipping")
    return(NULL)
  }
  tmp_sce = tmp_sce[, keep_ids]

  colData(tmp_sce)$zoom       = smpl_labels_dt[keep_ids]$zoom
  colData(tmp_sce)$sample_id  = sel_sample

  colData(tmp_sce)$cluster = factor(colData(tmp_sce)$zoom, levels = zoom_lbls)
  pb_label = agg_fn(tmp_sce, by_vars = c("cluster", "sample_id"),
    fun = "sum", all_cls = zoom_lbls)

  colData(tmp_sce)$cluster = factor("all_cells", levels = "all_cells")
  pb_all   = agg_fn(tmp_sce, by_vars = c("cluster", "sample_id"),
    fun = "sum", all_cls = "all_cells")

  missing_labels = setdiff(zoom_lbls, assayNames(pb_label))
  for (lbl in missing_labels) {
    assay(pb_label, lbl) = Matrix(0, nrow = nrow(pb_label), ncol = 1, sparse = FALSE,
      dimnames = list(rownames(pb_label), sel_sample))
  }

  list(
    assays  = c(as.list(assays(pb_label)), as.list(assays(pb_all))),
    rowData = rowData(tmp_sce)
  )
}


make_pseudobulk_celltypes <- function(config_f, out_dir = "output/compare_hvg_methods",
  label_col = "predicted_label_naive", min_prob = 0.8, n_cores = 8) {
  
  browser()
  p = .parse_config(config_f)
  
  # get zoom names and corresponding labels
  zoom_names_dt = .get_zoom_names_and_lbls(p$zoom) %>%
    setnames(old = 'label', new = 'predicted_label_naive')

  h5ads_f = file.path(p$proj_dir, "output", paste0(p$short_tag, "_integration"),
    sprintf("h5ads_clean_paths_%s.yaml", p$tag))
  h5ad_paths = read_yaml(h5ads_f)

  labels_f = file.path(p$proj_dir, "output", paste0(p$short_tag, "_label_celltypes"),
    sprintf("labels_scprocess_model_human_cns_%s.csv.gz", p$tag))
  
  # only keep cells that have confidently assigned zoom labels
  labels_dt = fread(labels_f) %>%
    .[probability_naive >= min_prob] %>%
    merge(zoom_names_dt, by = 'predicted_label_naive', all.x = TRUE) %>%
    .[!is.na(zoom)]
      
  sample_names = names(h5ad_paths)
  pb_ls = sample_names %>% lapply(FUN = .make_one_pseudobulk_celltypes,
    h5ad_paths = h5ad_paths, labels_dt = labels_dt,
    agg_fn = aggregateData_datatable)
  
  names(pb_ls) = sample_names
  pb_ls = Filter(Negate(is.null), pb_ls)

  assay_names_all = c(unique(zoom_names_dt$zoom), "all_cells")
  assay_ls = lapply(assay_names_all, function(nm) {
    mats = lapply(pb_ls, function(x) x$assays[[nm]])
    common_genes = rownames(mats[[1]])
    mats = lapply(mats, function(m) m[common_genes, , drop = FALSE])
    Reduce(cbind, mats)
  }) %>% setNames(assay_names_all)

  message("  merging pseudobulk counts")
  pb = SingleCellExperiment(assays = assay_ls, rowData = pb_ls[[1]]$rowData)
  colData(pb)$sample_id = colnames(pb)
 
  out_f = file.path(out_dir, p$proj_name, "pb_celltypes.rds")
  dir.create(dirname(out_f), recursive = TRUE, showWarnings = FALSE)
  message("  saving to ", out_f)
  saveRDS(pb, out_f, compress = FALSE)

  message("done!")
  return(NULL)
}
