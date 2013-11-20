#include "MultiCA_cpeigen.h"

using Rcpp::as;
using Rcpp::List;
using Rcpp::Named;
using Rcpp::wrap;
using Eigen::Lower;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::ArrayXd;
using Eigen::SelfAdjointEigenSolver; 
typedef Map<MatrixXd> MapMatd;
typedef Map<VectorXd> MapVecd;

MatrixXd AtA(const MapMatd& A) {
  int n(A.cols());
  return MatrixXd(n,n).setZero().selfadjointView<Lower>().rankUpdate(A.adjoint());
}

MatrixXd AAt(const MapMatd& A) {
  int n(A.rows());
  return MatrixXd(n,n).setZero().selfadjointView<Lower>().rankUpdate(A);
}


SEXP crossprodeigen(const SEXP M){
  const SelfAdjointEigenSolver<MatrixXd> es(AtA(as<MapMatd>(M)));  
  MatrixXd eivec(es.eigenvectors());
  ArrayXd eival(es.eigenvalues());
  return List::create(Named("eigenvectors") = eivec,
                      Named("eigenvalues") = eival);
}

SEXP tcrossprodeigen(const SEXP M){
  const SelfAdjointEigenSolver<MatrixXd> es(AAt(as<MapMatd>(M)));  
  MatrixXd eivec(es.eigenvectors());
  ArrayXd eival(es.eigenvalues());
  return List::create(Named("eigenvectors") = eivec,
                      Named("eigenvalues") = eival);
}



