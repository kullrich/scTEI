#include <RcppArmadillo.h>
#include <RcppThread.h>
#include <string.h>
#include <math.h>
#include <random>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' @useDynLib scTEI, .registration = TRUE
//' @import Rcpp
//' @import Matrix
//' @title rcpp_boottei_parallel
//' @name rcpp_boottei_parallel
//' @description computes the phylogenetically based
//' transcriptome evolutionary index (TEI) shuffling the strata for permutation
//' statistic
//' @return sparseMatrix
//' @param expression ExpressionSet as sparseMatrix
//' @param ps named Phylostratum
//' @param permutations number of permutations
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
//' ## get permutations
//' rcpp_boottei_parallel(spmat, ps, 100, 1)
//' @export rcpp_boottei_parallel
//' @author Kristian K Ullrich
NumericVector permut(const NumericVector& a) {
    NumericVector b = clone(a);
    std::random_device rng;
    std::mt19937 urng(rng());
    std::shuffle(b.begin(), b.end(), urng);
    return b;
}
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_boottei_parallel(const arma::sp_mat& expression,
    Rcpp::NumericVector ps,
    const int& permutations,
    int ncores = 1){
    std::vector< std::string > psnames =  ps.attr("names");
    int n_col = expression.n_cols;
    int n_row = expression.n_rows;
    Rcpp::NumericMatrix fMatrix(n_row, n_col);
    Rcpp::NumericVector sampledVector(ps.length());
    Rcpp::NumericMatrix sampledMatrix(ps.length(), permutations);
    for (size_t l = 0; l < permutations; l++) {
        sampledVector = permut(ps);
        for(size_t m = 0; m < ps.length(); m++) {
            sampledMatrix(m, l) = sampledVector[m];
        }
    }
    Rcpp::NumericVector sumx(n_col);
    Rcpp::NumericMatrix bootM(permutations,n_col);
    RcppThread::ProgressBar bar(n_col, 1);
    RcppThread::parallelFor(0, n_col, [&] (int j) {
        for (size_t i = 0; i < n_row; i++) {
            sumx[j] += expression(i, j);
        }
        for (size_t k = 0; k < n_row; k++) {
            fMatrix(k, j) = expression(k, j) / sumx[j];
        }
        for (size_t n = 0; n < permutations; n++) {
            double teisum = 0;
            for (size_t o = 0; o < n_row; o++) {
                teisum += (sampledMatrix(o, n) * fMatrix(o, j));
            }
            bootM(n, j) = teisum;
        }
    }, ncores);
    return bootM;
}
