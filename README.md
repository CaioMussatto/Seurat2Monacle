# Seurat2Monacle

**Seurat2Monacle**: Convert Seurat Object to Monocle CDS Seurat2Monacle is an R package that facilitates converting Seurat objects to the Monocle CellDataSet (CDS) format for single-cell analysis, including preprocessing, clustering, and UMAP visualization. The main function of the package, seurat2monacle(), allows you to transform a Seurat object (compatible with Seurat v3 and v5) into a CellDataSet from Monocle3.

------------------------------------------------------------------------

## Installation

### Install the package from GitHub

1.  **Install the necessary dependencies**:

To ensure that all required packages are installed, use the following code to install dependencies before installing the package:

``` r
# Install GitHub packages for monocle3 and SeuratWrappers 
install.packages("remotes") 
remotes::install_github("cole-trapnell-lab/monocle3") 
remotes::install_github("satijalab/seurat-wrappers") 
# Install the Seurat2Monacle package: Once the dependencies are installed, you can install your Seurat2Monacle package using devtools:
devtools::install_github("CaioMussatto/seurat2monacle")
```

2.  **Install the `Seurat2Monacle` package**:

Once the dependencies are installed, you can install your `Seurat2Monacle` package using `devtools`:

------------------------------------------------------------------------

## Usage

Load the library After installation, load the package into your R environment:

``` r
library(seurat2monacle)
```

### Main Function: **seurat2monacle()**

The `seurat2monacle()` function converts a `Seurat` object to a `CellDataSet` from `Monocle3`. You can specify a metadata column or feature for coloring the cells in the UMAP plot.

#### Parameters object:

-   `object`: A `Seurat` object to be converted.
-   `umap_color` (optional): The name of a metadata column or feature to color cells in the UMAP plot. If `NULL`, the default coloring will be used.

#### Example Usage

``` r
# Load Seurat package
library(Seurat)

# Assume you have a Seurat object called seurat_obj
# Here, we're using a fictitious example
seurat_obj <- readRDS("path/to/your_seurat_object.rds")

# Convert to Monocle CDS and generate a UMAP plot
cds <- seurat2monacle(seurat_obj, umap_color = "cell_type")

# The result is a CellDataSet object that can be used for further analysis
```

------------------------------------------------------------------------

#### Function Details

The `seurat2monacle()` function performs the following steps:

1.  **Checks if the assay in the Seurat object is in V5 format** and, if necessary, converts it to V3 format.

2.  **Extracts gene annotations, cell metadata, and the expression matrix** from the Seurat object.

3.  **Creates a new `CellDataSet` (CDS)** from the extracted data.

4.  **Assigns partitions and cluster labels** to the CDS based on Seurat's cell clusters.

5.  **Preprocesses and clusters the cells** in the CDS.

6.  **Generates a UMAP plot** with cells colored by the specified metadata column or feature.

------------------------------------------------------------------------

## Contributing

If you would like to contribute to the development of this package, feel free to submit a **pull request** or open an **issue** on GitHub.

------------------------------------------------------------------------

## License

This package is licensed under the MIT License.
