#ifndef _MultiCA_RCPP_CPEIGEN_H
#define _MultiCA_RCPP_CPEIGEN_H

#include <Rcpp.h>
#include <RcppEigen.h>

RcppExport SEXP crossprodeigen (const SEXP M) ;
RcppExport SEXP tcrossprodeigen (const SEXP M) ;

#endif