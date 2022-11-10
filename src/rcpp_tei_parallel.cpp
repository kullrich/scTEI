#include <RcppArmadillo.h>
#include <RcppThread.h>
#include <string.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' @useDynLib scTEI, .registration = TRUE
//' @import Rcpp
//' @import Matrix
//' @title rcpp_tei_parallel
//' @name rcpp_tei_parallel
//' @description computes the phylogenetically based
//' transcriptome evolutionary index (TEI)
//' @return list
//' @param expression ExpressionSet as sparseMatrix
//' @param ps named Phylostratum
//' @param ncores number of cores
//' @examples
//' ## load example PhyloExpressionSetExample
//'
//' data("PhyloExpressionSetExample", package="myTAI")
//'
//' ## convert into sparseMatrix - rownames GeneID
//'
//' spmat <- as(data.matrix(PhyloExpressionSetExample[,-c(1,2)]),
//'     "sparseMatrix")
//' rownames(spmat) <- PhyloExpressionSetExample$GeneID
//'
//' ## create named Phylostratum vector
//'
//' ps <- setNames(PhyloExpressionSetExample$Phylostratum,
//'     PhyloExpressionSetExample$GeneID)
//'
//' ## get TEI
//' rcpp_tei_parallel(spmat, ps)
//' @export rcpp_tei_parallel
//' @author Kristian K Ullrich
// [[Rcpp::export]]
Rcpp::List rcpp_tei_parallel(const arma::sp_mat& expression,
    Rcpp::NumericVector ps,
    int ncores = 1){
    std::vector< std::string > psnames =  ps.attr("names");
    int n_col = expression.n_cols;
    int n_row = expression.n_rows;
    Rcpp::NumericVector sumx(n_col);
    Rcpp::NumericVector teisum(n_col);
    Rcpp::NumericVector tei(n_col);
    RcppThread::ProgressBar bar(n_col, 1);
    RcppThread::parallelFor(0, n_col, [&] (int j) {
        for (size_t i = 0; i < n_row; i++) {
            sumx[j] += expression(i, j);
            teisum[j] += expression(i, j) * ps[i];
            tei[j] = teisum[j]/sumx[j];
        }
    }, ncores);
    return Rcpp::List::create(Rcpp::Named("sumx") = sumx,
        Rcpp::Named("teisum") = teisum, Rcpp::Named("tei") = tei);
}
