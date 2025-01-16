.onLoad <- function(libname, pkgname) {
  required_packages <- c(
    "BiocManager", "SeuratObject", "Seurat", "dplyr", "Matrix",
    "BiocGenerics", "DelayedArray", "DelayedMatrixStats", "limma", "lme4", "S4Vectors",
    "SingleCellExperiment", "SummarizedExperiment", "batchelor", "HDF5Array", "terra", "ggrastr"
  )

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  bioc_packages <- required_packages[required_packages %in% c(
    'BiocGenerics', 'DelayedArray', 'DelayedMatrixStats', 'limma', 'lme4', 'S4Vectors',
    'SingleCellExperiment', 'SummarizedExperiment', 'batchelor', 'HDF5Array', 'terra', 'ggrastr'
  )]

  if (length(bioc_packages) > 0) {
    BiocManager::install(bioc_packages)
  }

  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }

  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }

  remotes::install_github('cole-trapnell-lab/monocle3')

  remotes::install_github('satijalab/seurat-wrappers')


  invisible(lapply(required_packages, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Instalando e carregando", pkg))
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE, quietly = TRUE)
  }))
}
