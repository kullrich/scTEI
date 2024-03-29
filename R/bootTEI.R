#' @title Compute a Permutation Matrix of Transcriptome Evolutionary Index (TEI)
#' @description
#' This function computes the transcriptome evolutionary
#' index (TEI) using permuted strata values.
#' @param ExpressionSet expression object with rownames as GeneID (dgCMatrix)
#' or standard PhyloExpressionSet object.
#' @param Phylostratum a named vector representing phylostratum per GeneID with
#' names as GeneID (not used if Expression is PhyloExpressionSet).
#' @param permutations a numeric value specifying the number of permutations
#' to be performed.
#' @param split specify number of columns to split
#' @param showprogress boolean if progressbar should be shown
#' @param threads specify number of threads
#' @importFrom utils txtProgressBar
#' @importFrom myTAI is.ExpressionSet
#' @details The strata values are sampled and the global \code{\link{TEI}}
#' is calculated N times.
#' @return a numeric matrix storing the TEI values based on permuted
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
#' ## calculate TEI permutation values
#' bootTEI(
#'     ExpressionSet=celegans@assays$RNA@counts,
#'     Phylostratum=ps_vec
#' )
#'
#' @export bootTEI
#' @author Kristian K Ullrich

bootTEI <- function(ExpressionSet,
    Phylostratum=NULL,
    permutations=100,
    split=100000,
    showprogress=TRUE,
    threads=1){
    if(is(ExpressionSet, "Matrix")){
        common_ids<-sort(Reduce(intersect, list(rownames(ExpressionSet),
            names(Phylostratum))))
        Phylostratum<-Phylostratum[names(Phylostratum) %in% common_ids]
        Phylostratum<-Phylostratum[order(names(Phylostratum))]
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
            es_bm<-NULL
            for(i in seq_along(split_start)){
                es<-ExpressionSet[,split_start[i]:split_end[i]]
                es<-es[rownames(es) %in% common_ids,,drop=FALSE]
                es<-es[order(rownames(es)), ]
                es_bm<-rcpp_boottei_parallel(es, Phylostratum,
                    permutations, threads)
                colnames(es_bm)<-colnames(es)
                if(showprogress){
                    setTxtProgressBar(pb, i)
                }
                es <- NULL
                OUT<-cbind(OUT, es_bm)
            }
            bm<-OUT
        }else{
            ExpressionSet<-ExpressionSet[rownames(ExpressionSet) %in%
                common_ids,,drop=FALSE]
            ExpressionSet<-ExpressionSet[order(rownames(ExpressionSet)), ]
            bm<-rcpp_boottei_parallel(ExpressionSet, Phylostratum,
                permutations, threads)
            colnames(bm)<-colnames(ExpressionSet)
        }
    }
    if(is(ExpressionSet, "data.frame") | is(ExpressionSet, "tibble")){
        if(myTAI::is.ExpressionSet(ExpressionSet)){
            Phylostratum<-setNames(ExpressionSet$Phylostratum,
                ExpressionSet$GeneID)
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
            es_bm<-NULL
            for(i in seq_along(split_start)){
                es<-ExpressionSet[,split_start[i]:split_end[i]]
                es_bm<-rcpp_boottei_parallel(es, Phylostratum,
                    permutations, threads)
                colnames(es_bm)<-colnames(es)
                if(showprogress){
                    setTxtProgressBar(pb, i)
                }
                es<-NULL
                OUT<-cbind(OUT, es_bm)
            }
            bm<-OUT
        }else{
            bm<-rcpp_boottei_parallel(ExpressionSet, Phylostratum,
                permutations, threads)
            colnames(bm)<-colnames(ExpressionSet)
        }
    }
    return(bm)
}
