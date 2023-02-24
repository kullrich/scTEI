#' @title Transfer Seurat object to monocle3 object
#' @description
#' This function transfer Seurat object including PCA, tSNE, UMAP into monocole3
#' object.
#' @param seurat_obj Seurat object
#' @param seurat_assay Seurat assay
#' @param seurat_reduction Seurat reduction
#' @return cell data set object
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
#' 
#' ## convert into monocle3
#' celegans_cds<-seurat2monocle3(
#'     seurat_obj=celegans,
#'     seurat_assay="RNA",
#'     seurat_reduction=NULL
#' )
#' @export seurat2monocle3
#' @author Kristian K Ullrich

seurat2monocle3 <- function(seurat_obj, seurat_assay, seurat_reduction=NULL){
    Seurat::DefaultAssay(seurat_obj) <- seurat_assay
    cds.obj <- Seurat::as.SingleCellExperiment(x=seurat_obj, assay=seurat_assay)
    gene_info <- data.frame(gene_name=rownames(seurat_obj),
        gene_short_name=rownames(seurat_obj))
    rownames(gene_info) <- rownames(seurat_obj)
    SummarizedExperiment::rowData(cds.obj) <- gene_info
    cds <- as(object=cds.obj, Class='cell_data_set')
    SingleCellExperiment::reducedDimNames(cds)[
        which(SingleCellExperiment::reducedDimNames(cds)=="TSNE")] <- "tSNE"
    cds <- monocle3::estimate_size_factors(cds)
    if(!is.null(seurat_reduction)){
        if(seurat_reduction=="pca"){
            loadings <- Seurat::Loadings(seurat_obj@reductions$pca)
            stdev <- Seurat::Stdev(object=seurat_obj@reductions$pca)
        }
        if(seurat_reduction=="harmony"){
            loadings <- Seurat::Loadings(seurat_obj@reductions$harmony)
            stdev <- Seurat::Stdev(object=seurat_obj@reductions$harmony)
        }
        slot(object=cds, name='preprocess_aux')[['gene_loadings']] <- loadings
        slot(object=cds, name='preprocess_aux')[['prop_var_expl']] <- stdev
    }
    return(cds)
}
