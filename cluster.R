jhHeatmap <- function (dat_mat, pdf_fn, width = 12, height = 6, Colv = TRUE, 
          dendrogram = "both") {
  dat_mat <- as.matrix(dat_mat)
  if (!Colv) {
    dendrogram = "row"
  }
  pdf(pdf_fn, width = width, height = height)
  hm <- heatmap.2(dat_mat, Colv = Colv, dendrogram = dendrogram, 
                  col = heatmap_color("bluered", 5), hclust = function(x) hclust(x, 
                                                                                 method = "ward.D"), distfun = function(x) as.dist((1 - 
                                                                                                                                      cor(t(x)))/2), scale = "row", key = TRUE, symkey = FALSE, 
                  density.info = "none", trace = "none", cexRow = 0.5)
  dev.off()
  return(hm)
}
quick_heatmap_eda <- function (dat, top_var_percent = 0.9, outdir = "tall", ha = NULL, 
          heatmap_title = "leukemia", use_raster = F, column_split = 10, 
          height = 6, width = 16, cluster_columns = TRUE, show_row_names = F, 
          show_column_names = F, fig_type = "pdf", out_put_rds = F, 
          main_color = circlize::colorRamp2(c(-1.8, 0, 1.8), c("blue", 
                                                               "white", "red")), hdf = NULL, keep_annotation = F, column_order = NULL, 
          marker_gene = "", ...) {
  dat_var <- apply(dat, 1, var)
  checkdir(outdir)
  message(glue("We are processing percent: {top_var_percent}"))
  index <- dat_var >= quantile(dat_var, probs = top_var_percent, 
                               na.rm = T)
  if (top_var_percent == 0) {
    exp_plot <- data.matrix(dat)
  }
  else {
    exp_plot <- data.matrix(dat[index, ])
  }
  marker_gene <- base::intersect(marker_gene, rownames(exp_plot))
  if (is.null(marker_gene) || length(marker_gene) < 1) {
    rha <- NULL
  }
  else {
    atv <- lapply(marker_gene, function(gene_i) {
      which(rownames(exp_plot) == gene_i)
    }) %>% unlist()
    rha <- rowAnnotation(marker_gene = anno_mark(at = atv, 
                                                 labels = marker_gene))
  }
  hm <- jhHeatmap(as.matrix(exp_plot), glue("{outdir}/jhHeatmap_{top_var_percent}.pdf"), 
                  width = width, height = height, Colv = T)
  mat_scaled <- t(hm$carpet)[row.names(exp_plot), colnames(exp_plot)]
  rowdendrogram <- hm$rowDendrogram
  if (cluster_columns) {
    ht <- ComplexHeatmap::Heatmap(mat_scaled, name = heatmap_title, 
                                  use_raster = use_raster, col = main_color, row_order = rev(hm$rowInd), 
                                  column_order = hm$colInd, cluster_rows = rev(hm$rowDendrogram), 
                                  cluster_columns = hm$colDendrogram, show_row_names = show_row_names, 
                                  show_column_names = show_column_names, top_annotation = ha, 
                                  column_split = column_split, column_gap = unit(0.3, 
                                                                                 "mm"), column_title = NULL, right_annotation = rha, 
                                  ...)
    coldendrogram <- hm$colDendrogram
  }
  else {
    if (is.null(column_order)) {
      column_order <- 1:ncol(mat_scaled)
    }
    message("We cannot split column ...")
    ht <- ComplexHeatmap::Heatmap(mat_scaled, name = heatmap_title, 
                                  use_raster = use_raster, col = main_color, row_order = rev(hm$rowInd), 
                                  cluster_rows = rev(hm$rowDendrogram), cluster_columns = F, 
                                  clustering_method_columns = "ward.D2", show_row_names = show_row_names, 
                                  show_column_names = show_column_names, top_annotation = ha, 
                                  column_gap = unit(0.3, "mm"), column_title = NULL, 
                                  column_order = column_order, right_annotation = rha, 
                                  ...)
  }
  fgn <- glue("{outdir}/heatmap_{top_var_percent}.{fig_type}")
  create_fig <- glue("{fig_type}(fgn, height = {height}, width = {width})")
  eval(parse(text = create_fig))
  ht <- ComplexHeatmap::draw(ht, annotation_legend_side = "left", 
                             heatmap_legend_side = "right")
  dev.off()
  if (!cluster_columns) {
    message("No cluster column coldendrogram is NA ...")
    coldendrogram <- NA
  }
  cl <- map_int(column_order(ht), length)
  ci <- NULL
  for (i in 1:length(cl)) {
    ci <- c(ci, paste0(rep("G", cl[i]), i))
  }
  output_info <- data.frame(sample_id = colnames(mat_scaled)[unlist(column_order(ht))], 
                            sub_groups = ci)
  if (!is.null(hdf) & keep_annotation) {
    output_info <- cbind(output_info, hdf[unlist(column_order(ht)), 
    ])
  }
  heatmapeGenes <- data.frame(gene_name = rownames(mat_scaled)[unlist(row_order(ht))])
  if (out_put_rds) {
    rds_list <- list(complexheatmap = ht, heatmap2 = hm, 
                     rowdendrogram = rowdendrogram, coldendrogram = coldendrogram)
    write_rds(rds_list, file = glue("{outdir}/heatmap_{top_var_percent}_groups{column_split}.rds"))
  }
  list2excel(list(sampleinfo = output_info, heatmapeGenes = heatmapeGenes), 
             glue("{outdir}/sampleinfo_{top_var_percent}.xlsx"))
}