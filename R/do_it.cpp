#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List subsampling_nm_cpp(const arma::mat& X, const arma::vec& Y, const arma::mat& MN_plt, const arma::vec& P_plt, const std::string& criterion) {
  arma::vec nm;
  arma::vec norm;
  if (criterion == "OptA") {
    arma::mat invMN_X = arma::solve(MN_plt, X.t()).t();
    norm = arma::sqrt(arma::sum(invMN_X % invMN_X, 1));
    nm = arma::abs(Y - P_plt) % norm;
  } else if (criterion == "OptL") {
    norm = arma::sqrt(arma::sum(X % X, 1));
    nm = arma::abs(Y - P_plt) % norm;
  } else if (criterion == "LCC") {
    nm = arma::abs(Y - P_plt);
  }
  return Rcpp::List::create(Rcpp::Named("norm") = norm, Rcpp::Named("nm") = nm);
}
