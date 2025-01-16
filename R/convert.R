#' Convert Seurat Object to Monocle CDS
#'
#' This function converts a Seurat object into a Monocle CellDataSet (CDS).
#' It supports both Seurat V3 and V5 assays, automatically converting V5 assays
#' to the V3 format before processing. The resulting CDS is preprocessed and
#' clustered, and a UMAP plot is generated.
#'
#' @param object A Seurat object to be converted.
#' @param umap_color (Optional) The name of a metadata column or feature to
#'   use for coloring cells in the UMAP plot. If \code{NULL}, the default
#'   coloring will be used.
#'
#' @return A Monocle \code{CellDataSet} object containing the processed data.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Checks if the assay in the Seurat object is in V5 format and converts it to V3 if necessary.
#'   \item Extracts gene annotations, cell metadata, and the expression matrix from the Seurat object.
#'   \item Initializes a Monocle \code{CellDataSet} with the extracted data.
#'   \item Clusters cells and assigns partitions and cluster labels to the CDS.
#'   \item Embeds UMAP coordinates into the CDS.
#'   \item Optionally plots the UMAP with cells colored by the specified feature or metadata.
#' }
#'
#' @examples
#' \dontrun{
#' # Load Seurat object
#' seurat_obj <- readRDS("path/to/seurat_object.rds")
#'
#' # Convert to Monocle CDS
#' cds <- seurat2monacle(seurat_obj, umap_color = "cell_type")
#' }
#'
#' @import SeuratObject
#' @importFrom Matrix Matrix
#' @importFrom monocle3 new_cell_data_set preprocess_cds cluster_cells plot_cells
#' @export
seurat2monacle <- function(object, umap_color = NULL) {
  if (!inherits(object, "Seurat")) {
    stop("Input must be a Seurat object.")
  }

  if (!"pca" %in% names(object@reductions)) {
    stop("The Seurat object must contain PCA reduction.")
  }
  if (!"umap" %in% names(object@reductions)) {
    stop("The Seurat object must contain UMAP reduction.")
  }

  if ("Assay5" %in% class(object[["RNA"]])) {
    message("Converting Assay5 to Assay (Seurat V3).")
    object[["RNA"]] <- as(object = object[["RNA"]], Class = "Assay")

  }

  object[["RNA"]] <- as(object = object[["RNA"]], Class = "Assay")

  gene_annotation <- data.frame(gene_short_name = rownames(object@reductions[["pca"]]@feature.loadings),
                                row.names = rownames(object@reductions[["pca"]]@feature.loadings))

  cell_metadata <- data.frame(barcode = object@assays[["RNA"]]@counts@Dimnames[[2]],
                              row.names = object@assays[["RNA"]]@counts@Dimnames[[2]])
  cell_metadata <- cbind(cell_metadata, object@meta.data)


  expression_matrix <- object@assays[["RNA"]]@counts[rownames(object@reductions[["pca"]]@feature.loadings), ]


  cds_from_seurat <- monocle3::new_cell_data_set(expression_matrix,
                                                 cell_metadata = cell_metadata,
                                                 gene_metadata = gene_annotation)


  recreate.partition <- factor(rep(1, length(cds_from_seurat@colData@rownames)),
                               levels = unique(rep(1, length(cds_from_seurat@colData@rownames))))
  names(recreate.partition) <- cds_from_seurat@colData@rownames
  cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

  list_cluster <- object@active.ident
  names(list_cluster) <- object@assays[["RNA"]]@data@Dimnames[[2]]
  cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
  cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"

  cds_from_seurat@int_colData@listData$reducedDims@listData[["UMAP"]] <- object@reductions[["umap"]]@cell.embeddings


  cds_from_seurat <- monocle3::preprocess_cds(cds_from_seurat, num_dim = 100)
  cds_from_seurat <- monocle3::cluster_cells(cds_from_seurat, resolution = 1e-5)


  if (is.null(umap_color)) {
  print(plot_cells(cds_from_seurat,
                   label_groups_by_cluster = TRUE,
                   label_branch_points = TRUE,
                   graph_label_size = 4))
} else {
  print(plot_cells(cds_from_seurat,
                   color_cells_by = umap_color,
                   label_groups_by_cluster = TRUE,
                   label_branch_points = TRUE,
                   graph_label_size = 4))
}

  return(cds_from_seurat)
}
