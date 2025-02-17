% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert.R
\name{seurat2monacle}
\alias{seurat2monacle}
\title{Convert Seurat Object to Monocle CDS}
\usage{
seurat2monacle(object, umap_color = NULL)
}
\arguments{
\item{object}{A Seurat object to be converted.}

\item{umap_color}{(Optional) The name of a metadata column or feature to
use for coloring cells in the UMAP plot. If \code{NULL}, the default
coloring will be used.}
}
\value{
A Monocle \code{CellDataSet} object containing the processed data.
}
\description{
This function converts a Seurat object into a Monocle CellDataSet (CDS).
It supports both Seurat V3 and V5 assays, automatically converting V5 assays
to the V3 format before processing. The resulting CDS is preprocessed and
clustered, and a UMAP plot is generated.
}
\details{
The function performs the following steps:
\itemize{
  \item Checks if the assay in the Seurat object is in V5 format and converts it to V3 if necessary.
  \item Extracts gene annotations, cell metadata, and the expression matrix from the Seurat object.
  \item Initializes a Monocle \code{CellDataSet} with the extracted data.
  \item Clusters cells and assigns partitions and cluster labels to the CDS.
  \item Embeds UMAP coordinates into the CDS.
  \item Optionally plots the UMAP with cells colored by the specified feature or metadata.
}
}
\examples{
\dontrun{
# Load Seurat object
seurat_obj <- readRDS("path/to/seurat_object.rds")

# Convert to Monocle CDS
cds <- seurat2monacle(seurat_obj, umap_color = "cell_type")
}

}
