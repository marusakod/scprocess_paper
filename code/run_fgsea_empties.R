#!/usr/bin/env Rscript
#
# Run fgsea (GO BP) on edgeR empty-vs-cells results for each lineage
# (all cells + zoom subsets). Reuses .calc_enrichment() from scprocess.

suppressPackageStartupMessages({
  library(data.table)
  library(magrittr)
  library(stringr)
  library(assertthat)
  library(fgsea)
  library(yaml)
  library(BiocParallel)
  RhpcBLASctl::omp_set_num_threads(1L)
})

source("code/fgsea.R")

GSEA_CUT = 0.1

args = commandArgs(trailingOnly = TRUE)
assert_that(length(args) >= 4,
  msg = "Usage: Rscript run_fgsea_empties.R <config.yaml> <out_dir> <gsea_dir> <path_set> [zoom1 zoom2 ...]")

config_f = args[1]
out_dir  = args[2]
gsea_dir = args[3]
path_set = args[4]
zooms    = if (length(args) > 4) args[5:length(args)] else character(0)

assert_that(path_set %in% c("go_bp", "go_cc"),
  msg = "path_set must be 'go_bp' or 'go_cc'")

cfg        = read_yaml(config_f)
proj_dir   = cfg$project$proj_dir
short_tag  = cfg$project$short_tag
full_tag   = cfg$project$full_tag
date_stamp = cfg$project$date_stamp
ref_txome  = cfg$project$ref_txome
tag        = paste0(full_tag, "_", date_stamp)
proj_name  = basename(proj_dir)

species_prefix = if (ref_txome %in% c("human_2024", "human_2020")) {
  "c5"
} else if (ref_txome %in% c("mouse_2024", "mouse_2020")) {
  "m5"
} else {
  stop("Unsupported ref_txome: ", ref_txome)
}
species_suffix = if (species_prefix == "c5") "Hs" else "Mm"

gmt_f = file.path(gsea_dir, sprintf("%s.go.%s.v2023.1.%s.symbols.gmt",
  species_prefix, sub("go_", "", path_set), species_suffix))
assert_that(file.exists(gmt_f), msg = sprintf("GMT file not found: %s", gmt_f))
pathways = gmtPathways(gmt_f)


run_one_lineage <- function(edger_f, lineage) {
  message(sprintf("[%s] %s: %s", proj_name, lineage, edger_f))
  assert_that(file.exists(edger_f))

  edger_dt = fread(edger_f)
  edger_dt[, symbol := str_extract(gene_id, "^[^_]+")]

  min_fdr  = edger_dt[FDR > 0]$FDR %>% min
  ranked_dt = edger_dt[, .(
    symbol,
    gsea_val  = ifelse(FDR == 0, (log10(min_fdr) - 1) * sign(logFC) * -1,
                       log10(FDR) * sign(logFC) * -1),
    decreases = TRUE
  )]

  fgsea_dt = .calc_enrichment(ranked_dt, pathways, GSEA_CUT)
  if (is.null(fgsea_dt) || nrow(fgsea_dt) == 0) {
    message("  no results")
    return(NULL)
  }

  fgsea_dt = fgsea_dt[!is.na(pval)][, lineage := lineage][, path_set := path_set]
  return(fgsea_dt)
}


results = list()

# main pipeline (all cells)
main_edger_f = file.path(proj_dir, "output", paste0(short_tag, "_empties"),
  sprintf("edger_empty_genes_all_%s.csv.gz", tag))
res = run_one_lineage(main_edger_f, "all cells")
if (!is.null(res)) results = c(results, list(res))

# zoom lineages
for (z in zooms) {
  zoom_edger_f = file.path(proj_dir, "output", paste0(short_tag, "_zoom"), z,
    sprintf("edger_empty_genes_%s_%s_%s.csv.gz", full_tag, z, date_stamp))
  if (!file.exists(zoom_edger_f)) {
    message(sprintf("  [SKIP] zoom %s: file not found", z))
    next
  }
  res = run_one_lineage(zoom_edger_f, z)
  if (!is.null(res)) results = c(results, list(res))
}

if (length(results) == 0) {
  message("[DONE] no results to save")
  quit(status = 0)
}

combined_dt = rbindlist(results, fill = TRUE)

out_f = file.path(out_dir, proj_name, sprintf("fgsea_empties_%s.csv.gz", path_set))
dir.create(dirname(out_f), recursive = TRUE, showWarnings = FALSE)
fwrite(combined_dt, file = out_f)
message(sprintf("[DONE] saved %d rows to %s", nrow(combined_dt), out_f))
