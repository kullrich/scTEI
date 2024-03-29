% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootTEI.R
\name{bootTEI}
\alias{bootTEI}
\title{Compute a Permutation Matrix of Transcriptome Evolutionary Index (TEI)}
\usage{
bootTEI(
  ExpressionSet,
  Phylostratum = NULL,
  permutations = 100,
  split = 1e+05,
  showprogress = TRUE,
  threads = 1
)
}
\arguments{
\item{ExpressionSet}{expression object with rownames as GeneID (dgCMatrix)
or standard PhyloExpressionSet object.}

\item{Phylostratum}{a named vector representing phylostratum per GeneID with
names as GeneID (not used if Expression is PhyloExpressionSet).}

\item{permutations}{a numeric value specifying the number of permutations
to be performed.}

\item{split}{specify number of columns to split}

\item{showprogress}{boolean if progressbar should be shown}

\item{threads}{specify number of threads}
}
\value{
a numeric matrix storing the TEI values based on permuted
strata.
}
\description{
This function computes the transcriptome evolutionary
index (TEI) using permuted strata values.
}
\details{
The strata values are sampled and the global \code{\link{TEI}}
is calculated N times.
}
\examples{

## get Seurat object
celegans<-readRDS(file=system.file("extdata",
    "celegans.embryo.SeuratData.rds", package="scTEI")
)

## load Caenorhabditis elegans gene age estimation
celegans_ps<-readr::read_tsv(
   file=system.file("extdata",
   "Sun2021_Orthomap.tsv", package="scTEI")
)

## define Phylostratum
ps_vec<-setNames(
    as.numeric(celegans_ps$Phylostratum),
    celegans_ps$GeneID
)

## calculate TEI permutation values
bootTEI(
    ExpressionSet=celegans@assays$RNA@counts,
    Phylostratum=ps_vec
)

}
\references{
Domazet-Loso T. and Tautz D. (2010).
\emph{A phylogenetically based transcriptome age index mirrors ontogenetic
divergence patterns}. Nature (468): 815-818.

Quint M et al. (2012).
\emph{A transcriptomic hourglass in plant embryogenesis}.
Nature (490): 98-101.

Drost HG et al. (2015)
Mol Biol Evol. 32 (5): 1221-1231 doi:10.1093/molbev/msv012
}
\author{
Kristian K Ullrich
}
