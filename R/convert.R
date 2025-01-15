monocle_seurat <- function(object, umap_color = NULL){
  library(monocle3)
  library(Seurat)
  library(dplyr)
  library(SeuratWrappers)
  gene_annotation <- as.data.frame(rownames(object@reductions[["pca"]]@feature.loadings),
                                   row.names = rownames(object@reductions[["pca"]]@feature.loadings))

  colnames(gene_annotation) <- "gene_short_name"


  cell_metadata <- as.data.frame(object@assays[["RNA"]]@counts@Dimnames[[2]],
                                 row.names = object@assays[["RNA"]]@counts@Dimnames[[2]])
  colnames(cell_metadata) <- "barcode"

  cell_metadata <- cbind(cell_metadata, object@meta.data)


  New_matrix <- object@assays[["RNA"]]@counts
  New_matrix <- New_matrix[rownames(object@reductions[["pca"]]@feature.loadings), ]
  expression_matrix <- New_matrix

  cds_from_seurat <- new_cell_data_set(expression_matrix,
                                       cell_metadata = cell_metadata,
                                       gene_metadata = gene_annotation)


  recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
  names(recreate.partition) <- cds_from_seurat@colData@rownames
  recreate.partition <- as.factor(recreate.partition)

  cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition



  list_cluster <- object@active.ident
  names(list_cluster) <- object@assays[["RNA"]]@data@Dimnames[[2]]

  cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster



  cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"


  cds_from_seurat@int_colData@listData$reducedDims@listData[["UMAP"]] <-object@reductions[["umap"]]@cell.embeddings

  cds_from_seurat <- preprocess_cds(cds_from_seurat, num_dim = 100)

  cds_from_seurat <- cluster_cells(cds_from_seurat, resolution=1e-5)
  print(plot_cells(cds_from_seurat,
                   color_cells_by = umap_color,
                   label_groups_by_cluster=TRUE,
                   label_branch_points=TRUE,
                   graph_label_size=4))
  gc()
  return(cds_from_seurat)
}
