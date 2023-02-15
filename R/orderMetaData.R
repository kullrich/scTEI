#' @title Order meta.data from Seurat object according to data 
#' @description
#' This function order the meta.data from Seurat object according to Seurat
#' data.
#' @param seurat_obj Seurat object
#' @param seurat_data Seurat data
#' @return Seurat object
#' @examples
#'
#' ## get Seurat object
#' celegans<-readRDS(file=system.file("extdata",
#'     "celegans.embryo.SeuratData.rds", package="scTEI")
#' )
#' 
#' ## re-order meta.data according to cell order
#' celegans<-orderMetaData(
#'     seurat_obj=celegans,
#'     seurat_data=celegans@assays$RNA@data
#' )
#' @export orderMetaData
#' @author Kristian K Ullrich

orderMetaData <- function(seurat_obj, seurat_data){
    seurat_obj@meta.data <- seurat_obj@meta.data[
    match(colnames(seurat_data),
    rownames(seurat_obj@meta.data)), ]
    return(seurat_obj)
}
