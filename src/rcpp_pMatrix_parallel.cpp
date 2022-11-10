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
//' @title rcpp_pMatrix_parallel
//' @name rcpp_pMatrix_parallel
//' @description computes the partial
//' transcriptome evolutionary index (TEI) values for each single gene
//' @return sparseMatrix
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
//' ## get pMatrix
//' rcpp_pMatrix_parallel(spmat, ps)
//' @export rcpp_pMatrix_parallel
//' @author Kristian K Ullrich
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_pMatrix_parallel(const arma::sp_mat& expression,
    Rcpp::NumericVector ps,
    int ncores = 1){
    std::vector< std::string > psnames =  ps.attr("names");
    int n_col = expression.n_cols;
    int n_row = expression.n_rows;
    Rcpp::NumericMatrix pMatrix(n_row, n_col);
    Rcpp::NumericVector sumx(n_col);
    RcppThread::ProgressBar bar(n_col, 1);
    RcppThread::parallelFor(0, n_col, [&] (int j) {
        for (size_t i = 0; i < n_row; i++) {
            sumx[j] += expression(i, j);
        }
        for (size_t k = 0; k < n_row; k++) {
            pMatrix(k, j) = expression(k, j) * ps[k] / sumx[j];
        }
    }, ncores);
    return pMatrix;
}
