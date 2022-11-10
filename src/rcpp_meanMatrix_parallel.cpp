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
//' @title rcpp_meanMatrix_parallel
//' @name rcpp_meanMatrix_parallel
//' @description computes mean expression
//' @return sparseMatrix
//' @param expression ExpressionSet as sparseMatrix
//' @param ps named Phylostratum
//' @param psgroup ordered unique Phylostratum
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
//' psgroup <- sort(unique(ps))
//'
//' ## get meanMatrix
//' rcpp_meanMatrix_parallel(spmat, ps, psgroup)
//' @export rcpp_meanMatrix_parallel
//' @author Kristian K Ullrich
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_meanMatrix_parallel(const arma::sp_mat& expression,
    Rcpp::NumericVector ps,
    Rcpp::NumericVector psgroup,
    int ncores = 1){
    std::vector< std::string > psnames =  ps.attr("names");
    int n_col = expression.n_cols;
    int n_row = expression.n_rows;
    int n_psgroup = psgroup.length();
    Rcpp::NumericMatrix meanMatrix(n_psgroup, n_col);
    Rcpp::NumericVector sumx(n_col);
    Rcpp::NumericVector eColumn(n_row);
    std::unordered_map<double, double> index_psgroup;
    std::unordered_map<double, double> count_psgroup;
    for (size_t p = 0; p < n_psgroup; p++) {
        index_psgroup[psgroup[p]] = p;
    }
    for(size_t i = 0; i < n_row; i++) {
        count_psgroup[index_psgroup[ps[i]]] += 1;
    }
    RcppThread::ProgressBar bar(n_col, 1);
    RcppThread::parallelFor(0, n_col, [&] (int j) {
        for (size_t k = 0; k < n_row; k++) {
            meanMatrix(index_psgroup[ps[k]], j) +=
                expression(k, j)/count_psgroup[index_psgroup[ps[k]]];
        }
    }, ncores);
    return meanMatrix;
}
