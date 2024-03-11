#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat dL_sq_test(const arma::mat& X, const arma::mat& Y_matrix,
                     const arma::mat& P, const arma::vec& p, int K, int d) {

  arma::mat S = Y_matrix - P;
  arma::vec p_squared = arma::square(p);
  arma::mat Psi(K*d, K*d, arma::fill::zeros);

  for (unsigned int i = 0; i < p.n_elem; ++i) {
    arma::mat kron_term = arma::kron(S.row(i).t() * S.row(i), X.row(i).t() * X.row(i));
    Psi += kron_term / p_squared(i);
  }

  return Psi;
}
