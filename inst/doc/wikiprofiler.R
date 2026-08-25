## -----------------------------------------------------------------------------
#| include: false
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 11,
  fig.height = 7,
  dpi = 144,
  warning = FALSE,
  message = FALSE
)


## -----------------------------------------------------------------------------
library(wikiprofiler)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(knitr)


## -----------------------------------------------------------------------------
data(geneList, package = "DOSE")

de <- names(geneList)[1:100]
wp_res <- enrichWP(de, organism = "Homo sapiens")
wp_tbl <- as.data.frame(wp_res)

de_symbol <- bitr(
  de,
  fromType = "ENTREZID",
  toType = "SYMBOL",
  OrgDb = "org.Hs.eg.db"
)

value <- stats::setNames(geneList[de_symbol[, 1]], de_symbol[, 2])
pathway_id <- wp_tbl$ID[1]

kable(head(wp_tbl[, c("ID", "Description", "p.adjust")], 5), digits = 4)


## -----------------------------------------------------------------------------
wpplot(pathway_id) |>
  wp_bgfill(
    value = value,
    low = "darkgreen",
    high = "firebrick",
    legend_x = 0.88,
    legend_y = 0.95
  ) |>
  wp_shadowtext(bg.r = 2, bg.col = "white")


## -----------------------------------------------------------------------------
p_single <- wpplot(pathway_id) |>
  wp_bgfill(
    value = value,
    low = "darkgreen",
    high = "firebrick",
    legend_x = 0.88,
    legend_y = 0.95
  ) |>
  wp_shadowtext()

single_png <- file.path(tempdir(), "wikiprofiler-single-demo.png")
wpsave(p_single, single_png, width = 11, height = 7)
single_png


## -----------------------------------------------------------------------------
expr_tbl <- de_symbol[1:40, c("ENTREZID", "SYMBOL")]
expr_tbl$score <- unname(geneList[expr_tbl$ENTREZID])

expr_tbl_dup <- expr_tbl[1:10, ]
expr_tbl_dup$score <- expr_tbl_dup$score * 0.5
expr_tbl2 <- rbind(expr_tbl, expr_tbl_dup)

mapped_value <- wp_map(
  expr_tbl2,
  value_col = "score",
  id_col = "ENTREZID",
  mapping = de_symbol[, c("ENTREZID", "SYMBOL")],
  mapping_from = "ENTREZID",
  mapping_to = "SYMBOL",
  aggregator = "mean"
)

head(mapped_value, 10)


## -----------------------------------------------------------------------------
mapping_table <- attr(mapped_value, "mapping_table")
kable(head(mapping_table, 10), digits = 4)


## -----------------------------------------------------------------------------
wpplot(pathway_id) |>
  wp_bgfill(
    value = mapped_value,
    low = "navy",
    high = "goldenrod",
    legend_x = 0.88,
    legend_y = 0.95
  ) |>
  wp_shadowtext()


## -----------------------------------------------------------------------------
control_value <- mapped_value
case_value <- mapped_value + rep(c(-0.6, 0.9), length.out = length(mapped_value))

p_compare <- wpplot(pathway_id) |>
  wp_comparefill(
    value = case_value,
    control = control_value,
    mode = "difference",
    low = "steelblue4",
    high = "darkorange2",
    legend_x = 0.88,
    legend_y = 0.95
  ) |>
  wp_shadowtext()

p_compare


## -----------------------------------------------------------------------------
kable(head(p_compare$comparison, 10), digits = 4)


## -----------------------------------------------------------------------------
#| eval: false
# wpplot(pathway_id) |>
#   wp_comparefill(
#     value = case_value,
#     control = control_value,
#     mode = "log2_ratio",
#     pseudocount = 1
#   ) |>
#   wp_shadowtext()


## -----------------------------------------------------------------------------
batch_dir <- file.path(tempdir(), "wikiprofiler-batch-demo")
if (dir.exists(batch_dir)) {
  unlink(batch_dir, recursive = TRUE)
}

batch_plots <- wp_render(
  pathway = wp_tbl[, c("ID", "Description")],
  value = mapped_value,
  n = 2,
  name_col = "Description",
  dir = batch_dir,
  file_ext = "png",
  filename_template = "{index}_{id}_{name}",
  shadowtext = TRUE,
  width = 11,
  height = 7
)

batch_files <- list.files(batch_dir, full.names = TRUE)
batch_files


## -----------------------------------------------------------------------------
names(batch_plots)


## -----------------------------------------------------------------------------
#| eval: false
# plots <- wp_render(
#   pathway = wp_res,
#   value = mapped_value,
#   n = 6,
#   shadowtext = TRUE
# )


## -----------------------------------------------------------------------------
#| eval: false
# wp_res <- enrichWP(gene_ids, organism = "Homo sapiens")
# 
# mapped_value <- wp_map(
#   expr_table,
#   value_col = "logFC",
#   id_col = "ENTREZID",
#   mapping = id_map,
#   mapping_from = "ENTREZID",
#   mapping_to = "SYMBOL",
#   aggregator = "mean"
# )
# 
# plots <- wp_render(
#   pathway = wp_res,
#   value = mapped_value,
#   n = 6,
#   dir = "wp_batch",
#   name_col = "Description",
#   filename_template = "{index}_{id}_{name}",
#   shadowtext = TRUE
# )


## -----------------------------------------------------------------------------
#| eval: false
# plots <- wp_render(
#   pathway = wp_res,
#   value = case_value,
#   control = control_value,
#   n = 6,
#   dir = "wp_compare_batch",
#   name_col = "Description",
#   filename_template = "{index}_{id}_{name}",
#   shadowtext = TRUE
# )

