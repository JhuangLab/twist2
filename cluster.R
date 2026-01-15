## The follow codes are for leukemia cluster. If you have any question, please contact with Jinyan Huang (hiekeen@gmail.com)
##Ensure you have the necessary packages installed:
# install.packages(c("gplots", "circlize", "glue", "fs", "readr", "purrr"))
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ComplexHeatmap")

library(gplots)
library(ComplexHeatmap)
library(circlize)
library(glue)
library(fs)
library(readr)
library(purrr)

#' Internal Helper: Generate a heatmap using gplots (heatmap.2)
#'
#' This function creates a classic heatmap using 1-Pearson correlation distance
#' and Ward linkage. It is primarily used here to calculate the clustering 
#' and scaling matrix for the main function.
#'
#' @param dat_mat Numeric matrix of expression data.
#' @param pdf_fn Character string. Output filename for the intermediate PDF.
#' @param width Numeric. Width of the PDF.
#' @param height Numeric. Height of the PDF.
#' @param Colv Logical. Whether to cluster columns.
#' @param dendrogram Character. "both", "row", "column", or "none".
#'
#' @return The heatmap object from gplots::heatmap.2 containing clustering info.
jhHeatmap <- function(dat_mat, pdf_fn, width = 12, height = 6, 
                      Colv = TRUE, dendrogram = "both") {
  
  dat_mat <- as.matrix(dat_mat)
  
  # Logic: If columns are not clustered, only show row dendrogram
  if (!Colv) {
    dendrogram <- "row"
  }
  
  # specific color palette (falling back to gplots::bluered if heatmap_color is missing)
  # You can revert to your custom function: col = heatmap_color("bluered", 5)
  my_col <- gplots::bluered(75) 
  
  # Open PDF device
  pdf(pdf_fn, width = width, height = height)
  
  # Generate heatmap using 1-correlation distance and ward.D clustering
  # Note: scale="row" z-scores the rows
  hm <- gplots::heatmap.2(
    dat_mat, 
    Colv = Colv, 
    dendrogram = dendrogram, 
    col = my_col, 
    hclust = function(x) hclust(x, method = "ward.D"), 
    distfun = function(x) as.dist((1 - cor(t(x))) / 2), 
    scale = "row", 
    key = TRUE, 
    symkey = FALSE, 
    density.info = "none", 
    trace = "none", 
    cexRow = 0.5
  )
  
  dev.off() # Close device
  
  return(hm)
}

#' Quick Exploratory Data Analysis (EDA) Heatmap
#'
#' This function filters the top variable genes, performs clustering using the
#' logic from `heatmap.2` (correlation distance), and visualizes the result
#' using `ComplexHeatmap`.
#'
#' @param dat Numeric matrix or data frame. Raw expression data.
#' @param top_var_percent Numeric (0-1). Percentile threshold for filtering highly variable genes.
#' @param outdir Character. Directory to save outputs.
#' @param ha ComplexHeatmap::HeatmapAnnotation object. Top column annotation.
#' @param heatmap_title Character. Title of the heatmap legend.
#' @param use_raster Logical. Whether to use rasterization for large matrices (speeds up rendering).
#' @param column_split Integer. Number of k-means groups or cuts for column splitting (if clustering is ON).
#' @param height Numeric. Figure height.
#' @param width Numeric. Figure width.
#' @param cluster_columns Logical. Whether to cluster columns.
#' @param show_row_names Logical. Show gene names?
#' @param show_column_names Logical. Show sample names?
#' @param fig_type Character. "pdf" or "png".
#' @param out_put_rds Logical. If TRUE, saves the heatmap object as an RDS file.
#' @param main_color Color mapping function (circlize::colorRamp2).
#' @param hdf Data frame. Metadata for columns (used for output table).
#' @param keep_annotation Logical. Whether to merge metadata into the output info table.
#' @param column_order Numeric/Character vector. Specific order for columns if cluster_columns is FALSE.
#' @param marker_gene Character vector. List of genes to highlight on the right side.
#' @param ... Additional arguments passed to ComplexHeatmap::Heatmap.
#'
#' @return Invisibly returns the ComplexHeatmap object.
#' @export
quick_heatmap_eda <- function(dat, 
                              top_var_percent = 0.9, 
                              outdir = "tall", 
                              ha = NULL, 
                              heatmap_title = "leukemia", 
                              use_raster = FALSE, 
                              column_split = 10, 
                              height = 6, 
                              width = 16, 
                              cluster_columns = TRUE, 
                              show_row_names = FALSE, 
                              show_column_names = FALSE, 
                              fig_type = "pdf", 
                              out_put_rds = FALSE, 
                              main_color = circlize::colorRamp2(c(-1.8, 0, 1.8), c("blue", "white", "red")), 
                              hdf = NULL, 
                              keep_annotation = FALSE, 
                              column_order = NULL, 
                              marker_gene = "", 
                              ...) {
  
  # 1. Create output directory if it doesn't exist
  fs::dir_create(outdir)
  message(glue::glue("Processing top variable genes: percentile {top_var_percent}"))
  
  # 2. Filter data based on variance
  dat_var <- apply(dat, 1, var, na.rm = TRUE)
  
  if (top_var_percent == 0) {
    # If 0, use all data
    exp_plot <- data.matrix(dat)
  } else {
    # Calculate quantile threshold
    cutoff <- quantile(dat_var, probs = top_var_percent, na.rm = TRUE)
    index <- dat_var >= cutoff
    exp_plot <- data.matrix(dat[index, , drop = FALSE])
  }
  
  # 3. Handle Marker Gene Annotation
  # Ensure marker genes actually exist in the filtered matrix
  valid_markers <- intersect(marker_gene, rownames(exp_plot))
  
  if (length(valid_markers) < 1) {
    rha <- NULL
  } else {
    # Find indices for marker genes
    atv <- which(rownames(exp_plot) %in% valid_markers)
    # Create Row Annotation (labels only for marker genes)
    rha <- ComplexHeatmap::rowAnnotation(
      marker_gene = ComplexHeatmap::anno_mark(
        at = atv, 
        labels = rownames(exp_plot)[atv]
      )
    )
  }
  
  # 4. Perform Clustering/Scaling via jhHeatmap (Helper)
  # This step generates a temporary PDF and calculates the Z-scores ('carpet') 
  # based on the correlation distance logic.
  hm_file <- glue::glue("{outdir}/jhHeatmap_{top_var_percent}.pdf")
  hm <- jhHeatmap(as.matrix(exp_plot), hm_file, width = width, height = height, Colv = TRUE)
  
  # Extract scaled matrix (z-scores) and dendrograms from the helper object
  # Note: hm$carpet is transposed in heatmap.2, so we transpose it back
  mat_scaled <- t(hm$carpet)[rownames(exp_plot), colnames(exp_plot)]
  rowdendrogram <- hm$rowDendrogram
  
  # 5. Prepare ComplexHeatmap Configuration
  if (cluster_columns) {
    coldendrogram <- hm$colDendrogram
    
    # Heatmap with column clustering
    ht <- ComplexHeatmap::Heatmap(
      mat_scaled,
      name = heatmap_title,
      use_raster = use_raster,
      col = main_color,
      # Sync row order with the helper heatmap
      row_order = rev(hm$rowInd), 
      cluster_rows = rev(hm$rowDendrogram),
      # Sync column order with the helper heatmap
      column_order = hm$colInd, 
      cluster_columns = hm$colDendrogram,
      show_row_names = show_row_names,
      show_column_names = show_column_names,
      top_annotation = ha,
      column_split = column_split,
      column_gap = unit(0.3, "mm"),
      column_title = NULL,
      right_annotation = rha,
      ...
    )
  } else {
    # Heatmap without column clustering
    message("Column clustering is OFF. Using provided order or default.")
    
    if (is.null(column_order)) {
      column_order <- 1:ncol(mat_scaled)
    }
    
    coldendrogram <- NA # No dendrogram
    
    ht <- ComplexHeatmap::Heatmap(
      mat_scaled,
      name = heatmap_title,
      use_raster = use_raster,
      col = main_color,
      row_order = rev(hm$rowInd),
      cluster_rows = rev(hm$rowDendrogram),
      cluster_columns = FALSE, # Explicitly FALSE
      show_row_names = show_row_names,
      show_column_names = show_column_names,
      top_annotation = ha,
      column_gap = unit(0.3, "mm"),
      column_title = NULL,
      column_order = column_order,
      right_annotation = rha,
      ...
    )
  }
  
  # 6. Save the Final Heatmap
  fgn <- glue::glue("{outdir}/heatmap_{top_var_percent}.{fig_type}")
  
  # Use match.fun instead of eval(parse) for safety
  # This looks for function "pdf" or "png" based on the string
  dev_fun <- match.fun(tolower(fig_type))
  
  dev_fun(fgn, height = height, width = width)
  
  # Draw the heatmap
  ht <- ComplexHeatmap::draw(ht, 
                             annotation_legend_side = "left", 
                             heatmap_legend_side = "right")
  dev.off()
  
  # 7. Extract Cluster Groups for Output
  # Get the column order after drawing
  final_col_order <- ComplexHeatmap::column_order(ht)
  
  # Handle cases where column_order might be a list (split) or vector
  if (is.list(final_col_order)) {
    cl <- purrr::map_int(final_col_order, length)
    # Create group labels (e.g., G1, G1, G2, G2...)
    ci <- rep(paste0("G", seq_along(cl)), times = cl)
    flat_col_indices <- unlist(final_col_order)
  } else {
    ci <- rep("G1", length(final_col_order))
    flat_col_indices <- final_col_order
  }
  
  output_info <- data.frame(
    sample_id = colnames(mat_scaled)[flat_col_indices],
    sub_groups = ci,
    stringsAsFactors = FALSE
  )
  
  # Merge with original metadata if requested
  if (!is.null(hdf) && keep_annotation) {
    # Ensure hdf rows align with sample IDs
    output_info <- cbind(output_info, hdf[flat_col_indices, , drop = FALSE])
  }
  
  # Extract ordered genes
  # row_order(ht) returns indices
  heatmapeGenes <- data.frame(
    gene_name = rownames(mat_scaled)[unlist(ComplexHeatmap::row_order(ht))]
  )
  
  # 8. Save Object to RDS (Optional)
  if (out_put_rds) {
    rds_list <- list(
      complexheatmap = ht,
      heatmap2 = hm,
      rowdendrogram = rowdendrogram,
      coldendrogram = coldendrogram,
      groups = output_info,
      genes = heatmapeGenes
    )
    rds_filename <- glue::glue("{outdir}/heatmap_{top_var_percent}_groups{column_split}.rds")
    readr::write_rds(rds_list, file = rds_filename)
    message(glue::glue("Saved RDS to: {rds_filename}"))
  }
  return(ht)
}