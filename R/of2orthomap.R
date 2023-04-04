#' @title Extract orthomap from OrthoFinder
#' @description
#' This function extract an orthomap from OrthoFinder results given a query
#' species.
#' @param seqname sequence name of the query species in OrthoFinder
#' @param qt query species taxid
#' @param sl species list as <orthofinder name><tab><species taxid>
#' @param oc specify OrthoFinder <Orthogroups.GeneCounts.tsv>
#' @param og specify OrthoFinder <Orthogroups.tsv>
#' @return orthomap object
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
#' @export of2orthomap
#' @author Kristian K Ullrich

of2orthomap <- function(seqname, qt, sl, oc, og){
    get_qtid <- function(qt){
        qt_name <- taxizedb::taxid2name(qt, db="ncbi")
        qlineage <- taxizedb::classification(qt, db="ncbi")
        qlineage_df <- data.frame(PSnum=rownames(qlineage[[1]]),
                              PStaxID=qlineage[[1]]$id,
                              PSname=qlineage[[1]]$name)
        qk <- NULL
        if(qlineage_df$PStaxID[2]=="2"){
            qk <- "Bacteria"
        }
        if(qlineage_df$PStaxID[2]=="2157"){
            qk <- "Archea"
        }
        if(qlineage_df$PStaxID[2]=="2759"){
            qk <- "Eukaryota"
        }
        out <- setNames(list(qt_name, qt, qlineage, qlineage_df, qk),
                        c("qt_name", "qt", "qlineage", "qlineage_df", "qk"))
        return(out)        
    }
    get_youngest_common <- function(ql, tl, get="id"){
        gyc <- tl[[1]][rev(which(tl[[1]]$id %in% ql[[1]]$id))[1],,drop=FALSE]
        if(get=="id"){
            return(gyc$id)
        }
        if(get=="name"){
            return(gyc$name)
        }
    }
    get_oldest_common <- function(ql, tl_ids){
        ql_ids <- ql[[1]]$id
        goc <- ql_ids[min(which(ql_ids %in% tl_ids))]
        return(goc)
    }
    qlineage <- get_qtid(qt)
    species_list <- readr::read_tsv(sl, col_names = c("species", "taxID"))
    species_lineages <- lapply(species_list$taxID, function(x) {
        taxizedb::classification(x, db="ncbi")})
    species_list$lineage <- species_lineages
    species_list$gyc <- unlist(lapply(species_list$lineage, function(x) {
        get_youngest_common(qlineage$qlineage, x)
    }))
    species_list$gyc_name <- unlist(lapply(species_list$lineage, function(x) {
        get_youngest_common(qlineage$qlineage, x, get="name")
    }))
    oc_df <- readr::read_tsv(oc, col_names=TRUE)
    oc_species <- colnames(oc_df)
    oc_seqname_idx <- which(colnames(oc_df)==seqname)
    oc_df <- oc_df[oc_df[, oc_seqname_idx]>0, , drop=FALSE]
    oc_og_hits <- apply(oc_df[, species_list$species], 1, function(x) {
        names(which(x!=0))
    })
    oc_df$goc <- unlist(lapply(oc_og_hits, function(x) {
        get_oldest_common(qlineage$qlineage,
                          species_list$gyc[which(species_list$species %in% x)])
    }))
    og_df <- readr::read_tsv(og, col_names=TRUE)
    og_species <- colnames(og_df)
    og_seqname_idx <- which(colnames(og_df)==seqname)
    og_df_query <- og_df[, c(1, og_seqname_idx)]
    og_df_query <- og_df_query[og_df_query$Orthogroup %in% oc_df$Orthogroup,
                               ,drop=FALSE]
    og_df_query$goc <- oc_df$goc
    omap <- data.frame(
        stringr::str_split_fixed(unlist(apply(og_df_query, 1, function(x) {
        paste(gsub(" ", "", unlist(strsplit(x[2], ","))), x[1],
              which(qlineage$qlineage[[1]]$id==x[3]), x[3])})),
        " ", 4))
    colnames(omap) <- c("seqID", "Orthogroup", "PSnum", "PStaxID")
    return(omap)
}
