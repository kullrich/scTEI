#' @title Get list of cell identity classes
#' @description
#' This function get a named list of cell identity classes.
#' @param cell_ident cell identity classes
#' @return a named list by identity class and corresponding cells
#' @examples
#'
#' ## get Seurat object
#' celegans<-readRDS(file=system.file("extdata",
#'     "celegans.embryo.SeuratData.rds", package="scTEI")
#' )
#' 
#' cell_groups<-Ident2cellList(Idents(celegans))
#' cell_groups
#' @export Ident2cellList
#' @author Kristian K Ullrich

Ident2cellList <- function(cell_ident){
    cL <- setNames(
    lapply(names(table(cell_ident)),
        function(x){which(cell_ident==x)}),
    names(table(cell_ident)))
    return(cL)
}
