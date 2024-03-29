% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seurat2monocle3.R
\name{seurat2monocle3}
\alias{seurat2monocle3}
\title{Transfer Seurat object to monocle3 object}
\usage{
seurat2monocle3(seurat_obj, seurat_assay, seurat_reduction = NULL)
}
\arguments{
\item{seurat_obj}{Seurat object}

\item{seurat_assay}{Seurat assay}

\item{seurat_reduction}{Seurat reduction}
}
\value{
cell data set object
}
\description{
This function transfer Seurat object including PCA, tSNE, UMAP into monocole3
object.
}
\examples{

## get Seurat object
celegans<-readRDS(file=system.file("extdata",
    "celegans.embryo.SeuratData.rds", package="scTEI")
)

## re-order meta.data according to cell order
celegans<-orderMetaData(
    seurat_obj=celegans,
    seurat_data=celegans@assays$RNA@data
)

## convert into monocle3
celegans_cds<-seurat2monocle3(
    seurat_obj=celegans,
    seurat_assay="RNA",
    seurat_reduction=NULL
)
}
\author{
Kristian K Ullrich
}
