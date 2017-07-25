#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector cara_cpp(NumericVector xs, double theta) {
  NumericVector rslt = xs; 
  for (int i = 0; i < xs.size(); i++) 
	rslt[i] = (1 - exp(-theta*xs[i]))/theta;
  return rslt;
}


// [[Rcpp::export]]
NumericVector w_prelec_cpp(NumericVector ps, double alpha, double beta=1) {
  NumericVector rslt = ps; 
  for (int i = 0; i < ps.size(); i++) 
	rslt[i] = exp(-beta * pow(-log(ps[i]), alpha));
  return rslt;
}


// [[Rcpp::export]]
NumericVector w_prelec_inv_cpp(NumericVector ps, double alpha, double beta=1) {
  NumericVector rslt = ps; 
  double gamma = 1/alpha;
  for (int i = 0; i < ps.size(); i++) 
	rslt[i] = exp(-pow(-log(ps[i])/beta, gamma));
  return rslt;
}


// [[Rcpp::export]]
NumericVector pow_seq(double delta, int n) {
  // return (delta, delta^2,..,delta^n)
  NumericVector rslt(n); 
  for (int i = 0; i < n; i++) 
	rslt[i] = pow(delta, i+1); 
  return rslt;
}
