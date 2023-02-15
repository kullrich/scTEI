#' @title Compute Relative Expression Levels
#' @description
#' This function computes the relative expression profiles.
#'
#' In detail, the \emph{Relative Expression Profile} is being computed as
#' follows over developmental stages s:
#'
#' \deqn{f_s = (e_s - e_min)/(e_max - e_min)},
#'
#' as follows over cells:
#'
#' \deqn{f_c = (e_c - e_min)/(e_max - e_min)},
#'
#' as follows over celltypes:
#'
#' \deqn{f_ct = (e_ct - e_min)/(e_max - e_min)},
#'
#' as follows over phylostrata:
#'
#' \deqn{f_ps = (e_ps - e_min)/(e_max - e_min)},
#'
#' where e_min and e_max denote either the minimum/maximum mean expression level
#' over developmental stages \eqn{s}, cells \eqn{c}, celltypes \eqn{ct} or
#' phylostrata \eqn{ps}.
#'
#' This linear transformation corresponds to a shift by \eqn{e_min - e_max}. As
#' a result, the relative expression level \eqn{f_s} of developmental stage
#' \eqn{s}, \eqn{f_c} of cell \eqn{c}, \eqn{f_ct} of celltype \eqn{ct} or
#' \eqn{f_ps} of phylotstratum \eqn{ps} with minimum \eqn{e_s}, \eqn{e_c},
#' \eqn{e_ct} or \eqn{e_ps} is 0, whereas the relative expression
#' level \eqn{f_s} of developmental stage \eqn{s}, \eqn{f_c} of cell \eqn{c},
#' \eqn{f_ct} of celltype \eqn{ct } or \eqn{f_ps} of phylotstratum \eqn{ps}
#' with maximum \eqn{e_s}, \eqn{e_c}, \eqn{e_ct} or \eqn{e_ps} is 1, and the
#' relative expression levels of all other stages \eqn{s}, cells \eqn{c},
#' celltypes \eqn{ct} or phylostrata \eqn{ps} range between 0 and 1.
#' @param ExpressionSet expression object with rownames as GeneID (dgCMatrix)
#' or standard PhyloExpressionSet object.
#' @param Phylostratum a named vector representing phylostratum per GeneID with
#' names as GeneID (not used if Expression is PhyloExpressionSet).
#' @param by specify min/max transformation by row (stages, cells, celltypes)
#' or by column (phylostratum)
#' @param groups specify stages or cells to be grouped into celltypes
#' by a named list
#' @param split specify number of columns to split
#' @param showprogress boolean if progressbar should be shown
#' @param threads specify number of threads
#' @importFrom utils txtProgressBar
#' @importFrom myTAI is.ExpressionSet
#' @details The partial TEI values combined per strata give an overall
#' impression of the contribution of each
#' strata to the global \code{\link{TEI}} pattern.
#' @return a numeric matrix storing the summed partial TEI values for each
#' strata.
#' @references
#' Domazet-Loso T. and Tautz D. (2010).
#' \emph{A phylogenetically based transcriptome age index mirrors ontogenetic
#' divergence patterns}. Nature (468): 815-818.
#'
#' Quint M et al. (2012).
#' \emph{A transcriptomic hourglass in plant embryogenesis}.
#' Nature (490): 98-101.
#'
#' Drost HG et al. (2015)
#' Mol Biol Evol. 32 (5): 1221-1231 doi:10.1093/molbev/msv012
#'
#' @examples
#'
#' ## get Seurat object
#' celegans<-readRDS(file=system.file("extdata",
#'     "celegans.embryo.SeuratData.rds", package="scTEI")
#' )
#'
#' ## load Caenorhabditis elegans gene age estimation
#' celegans_ps<-readr::read_tsv(
#'    file=system.file("extdata",
#'    "Sun2021_Orthomap.tsv", package="scTEI")
#' )
#'
#' ## define Phylostratum
#' ps_vec<-setNames(
#'     as.numeric(celegans_ps$Phylostratum),
#'     celegans_ps$GeneID
#' )
#'
#' ## get relative expression
#' Seurat::Idents(celegans)<-"embryo.time.bin"
#' reM<-REMatrix(
#'    ExpressionSet=celegans@assays$RNA@counts,
#'    Phylostratum=ps_vec
#' )
#'
#' ## get relative expression per cell group
#' Seurat::Idents(celegans)<-"embryo.time.bin"
#' cell_groups<-Ident2cellList(Idents(celegans))
#' reM<-REMatrix(
#'    ExpressionSet=celegans@assays$RNA@counts,
#'    Phylostratum=ps_vec,
#'    groups=cell_groups
#' )
#' p1<-ComplexHeatmap::Heatmap(
#'     reM,
#'     name="RExp",
#'     column_title="Relative Expression Profile - cell groups",
#'     row_title="More Recent <<< More Ancient",
#'     cluster_rows=FALSE,
#'     cluster_columns=FALSE,
#'     col=viridis::viridis(3)
#' )
#' p1
#' 
#' ## get relative expression over stages per cell group
#' Seurat::Idents(celegans)<-"embryo.time.bin"
#' cell_groups<-Ident2cellList(Idents(celegans))
#' reM<-REMatrix(
#'    ExpressionSet=celegans@assays$RNA@counts,
#'    Phylostratum=ps_vec,
#'    groups=cell_groups,
#'    by="row"
#' )
#' p2<-ComplexHeatmap::Heatmap(
#'     reM,
#'     name="RExp",
#'     column_title="Relative Expression Profile - cell groups - by row",
#'     row_title="More Recent <<< More Ancient",
#'     cluster_rows=FALSE,
#'     cluster_columns=FALSE,
#'     col=viridis::viridis(3)
#' )
#' p2
#'
#' ## get relative expression over phylostrata per cell group
#' Seurat::Idents(celegans)<-"embryo.time.bin"
#' cell_groups<-Ident2cellList(Idents(celegans))
#' reM<-REMatrix(
#'    ExpressionSet=celegans@assays$RNA@counts,
#'    Phylostratum=ps_vec,
#'    groups=cell_groups,
#'    by="column"
#' )
#' p3<-ComplexHeatmap::Heatmap(
#'     reM,
#'     name="RExp",
#'     column_title="Relative Expression Profile - cell groups - by col",
#'     row_title="More Recent <<< More Ancient",
#'     cluster_rows=FALSE,
#'     cluster_columns=FALSE,
#'     col=viridis::viridis(3)
#' )
#' p3
#' @export REMatrix
#' @author Kristian K Ullrich

REMatrix <- function(ExpressionSet,
    Phylostratum=NULL,
    by=NULL,
    groups=NULL,
    split=100000,
    showprogress=TRUE,
    threads=1){
    RE <- function(x){
        f_min<-min(x)
        f_max<-max(x)
        return((x-f_min)/(f_max-f_min))
    }
    if(is(ExpressionSet, "Matrix")){
        common_ids<-sort(Reduce(intersect, list(rownames(ExpressionSet),
            names(Phylostratum))))
        Phylostratum<-Phylostratum[names(Phylostratum) %in% common_ids]
        Phylostratum<-Phylostratum[order(names(Phylostratum))]
        PhylostratumGroups<-sort(unique(Phylostratum))
        if(ncol(ExpressionSet)>split){
            split_start<-seq(from=1, to=ncol(ExpressionSet), by=split)
            if(ncol(ExpressionSet)%%split!=0){
                split_end<-c(seq(from=split, to=ncol(ExpressionSet),
                    by=split), ncol(ExpressionSet))
            }else{
                split_end<-seq(from=split, to=ncol(ExpressionSet), by=split)
            }
            if(rev(split_end-split_start)[1]==0){
                split_start<-split_start[-length(split_start)]
                split_end[length(split_end)-1] <- rev(split_end)[1]
                split_end<-split_end[-length(split_end)]
            }
            if(showprogress){
                pb<-txtProgressBar(min=1, max=length(split_start), style=3)
            }
            OUT<-NULL
            es<-NULL
            es_meanMatrix<-NULL
            for(i in seq_along(split_start)){
                es<-ExpressionSet[,split_start[i]:split_end[i]]
                es<-es[rownames(es) %in% common_ids,,drop=FALSE]
                es<-es[order(rownames(es)), ]
                es<-ExpressionSet[,split_start[i]:split_end[i]]
                es_meanMatrix<-as(
                    rcpp_meanMatrix_parallel(es, Phylostratum,
                    PhylostratumGroups, threads), "sparseMatrix")
                colnames(es_meanMatrix)<-colnames(es)
                rownames(es_meanMatrix)<-PhylostratumGroups
                if(showprogress){
                    setTxtProgressBar(pb, i)
                }
                es<-NULL
                OUT<-cbind(OUT, es_meanMatrix)
            }
            meanMatrix<-OUT
        }else{
            ExpressionSet<-ExpressionSet[rownames(ExpressionSet) %in%
                common_ids,,drop=FALSE]
            ExpressionSet<-ExpressionSet[order(rownames(ExpressionSet)), ]
            meanMatrix<-as(
                rcpp_meanMatrix_parallel(ExpressionSet, Phylostratum,
                PhylostratumGroups, threads), "sparseMatrix")
            colnames(meanMatrix)<-colnames(ExpressionSet)
            rownames(meanMatrix)<-PhylostratumGroups
        }
    }
    if(is(ExpressionSet, "data.frame") | is(ExpressionSet, "tibble")){
        if(myTAI::is.ExpressionSet(ExpressionSet)){
            Phylostratum<-setNames(ExpressionSet$Phylostratum,
                ExpressionSet$GeneID)
            PhylostratumGroups<-sort(unique(Phylostratum))
            ExpressionSet<-as(data.matrix(
                ExpressionSet[,3:ncol(ExpressionSet)]), "sparseMatrix")
            rownames(ExpressionSet)<-names(Phylostratum)
        }
        if(ncol(ExpressionSet)>split){
            split_start<-seq(from=1, to=ncol(ExpressionSet), by=split)
            if(ncol(ExpressionSet)%%split != 0){
                split_end<-c(seq(from=split, to=ncol(ExpressionSet),
                    by=split), ncol(ExpressionSet))
            }else{
                split_end<-seq(from=split, to=ncol(ExpressionSet), by=split)
            }
            if(rev(split_end-split_start)[1]==0){
                split_start<-split_start[-length(split_start)]
                split_end[length(split_end)-1] <- rev(split_end)[1]
                split_end<-split_end[-length(split_end)]
            }
            if(showprogress){
                pb<-txtProgressBar(min=1, max=length(split_start), style=3)
            }
            OUT<-NULL
            es<-NULL
            es_meanMatrix<-NULL
            for(i in seq_along(split_start)){
                es<-ExpressionSet[,split_start[i]:split_end[i]]
                es_meanMatrix<-as(
                    rcpp_meanMatrix_parallel(es, Phylostratum,
                    PhylostratumGroups, threads), "sparseMatrix")
                colnames(es_meanMatrix)<-colnames(es)
                rownames(es_meanMatrix)<-PhylostratumGroups
                if(showprogress){
                    setTxtProgressBar(pb, i)
                }
                es<-NULL
                OUT<-cbind(OUT, es_meanMatrix)
            }
            meanMatrix<-OUT
        }else{
            meanMatrix<-as(
                rcpp_meanMatrix_parallel(ExpressionSet, Phylostratum,
                PhylostratumGroups, threads), "sparseMatrix")
            colnames(meanMatrix)<-colnames(ExpressionSet)
            rownames(meanMatrix)<-PhylostratumGroups
        }
    }
    if(is.null(groups)){
        if(!is.null(by)){
            if(by=="row"){
                meanMatrix<-t(apply(meanMatrix, 1, RE))
            }
            if(by=="column"){
                meanMatrix<-apply(meanMatrix, 2, RE)
            }
        }else{}
    }else{
        if(!is.null(by)){
            if(by=="row"){
                meanMatrix<-do.call(cbind, lapply(groups, function(x)
                    apply(meanMatrix[, x, drop=FALSE], 1, mean)))
                meanMatrix<-t(apply(meanMatrix, 1, RE))
            }
            if(by=="column"){
                meanMatrix<-do.call(cbind, lapply(groups, function(x)
                    apply(meanMatrix[, x, drop=FALSE], 1, mean)))
                meanMatrix<-apply(meanMatrix, 2, RE)
            }
        }else{
            meanMatrix<-do.call(cbind, lapply(groups, function(x)
                apply(meanMatrix[, x, drop=FALSE], 1, mean)))
        }
    }
    return(meanMatrix)
}
