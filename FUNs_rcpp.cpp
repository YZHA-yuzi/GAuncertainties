#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat comp_hazard(arma::mat h0){
  int nobs = h0.n_rows;
  int nweeks = h0.n_cols;
  arma::mat out;
  out.reshape(nobs, nweeks);
  for(int i=0; i<nobs; ++i){
    for(int j=0; j<nweeks; ++j){
      double prob_cond = 1;
      if(j == 0){
        prob_cond *= 1;
      }else{
        int kindex = j-1;
        for(int k=0; k<=kindex; ++k){
          prob_cond *= (1 - h0(i, k));
        }
      }
      out(i,j) = h0(i, j)*prob_cond;
    } //end loop over j; 27-36 weeks
  } // end loop over i; the number of observations
  return out;
}

