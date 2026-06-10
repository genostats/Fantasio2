#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// to optimize h(alpha) = sum_i log10( alpha * 10**flod[i] + 1 - alpha ) for alpha in [0,1]
// this is Newton method alpha = alpha - h'(alpha)/h''(alpha)
// noting Z = 10**flod[i] we have
// h'(alpha)  = 1/log(10) * sum_i (Z|i] - 1)/(alpha*Z[i] + 1 - alpha)
// h''(alpha) = -1/log(10) * sum_i [(Z|i] - 1)/(alpha*Z[i] + 1 - alpha)]^2
// h''(alpha) is always negative so the function is concave : only one max
// and if h'(0) < 0 the max is in 0, if h'(1) > 0 maximum is in 1
// [note: the 1/log(10) is dropped in the computations]
inline double sq(double x) {
  return x*x; 
}

// [[Rcpp::export]]
NumericVector hflod(NumericVector flod, double eps = 0.01) {
  int n = flod.size();
  std::vector<double> Z(n), Z1(n);
  for(int i = 0; i < n; i++) Z[i] = std::pow(10.0, flod[i]);
  // derivative in 0
  double d0 = 0;
  for(int i = 0; i < n; i++) d0 += Z[i] - 1;
  if(d0 < 0) // max for alpha = 0
    return NumericVector::create(0, 0);
  // derivative in 1
  double d1 = 0;
  for(int i = 0; i < n; i++) d1 += 1 - 1/Z[i];
  if(d1 > 0) {
    double hflod = 0;
    for(int i = 0; i < n; i++) hflod += flod[i];
    return NumericVector::create(1, hflod);
  }
  // Newton iterations, starting from 0
  double dd0 = 0; // in fact, the opposite of the 2nd derivative
  for(int i = 0; i < n; i++) dd0 += sq(Z[i] - 1);
  double alpha = d0/dd0;
  alpha = (alpha < 1)?alpha:1;
  double alpha1 = 0;
  while(true) {
    double da = 0;   // derivative in alpha
    double dda = 0;  // opposite of the second derivative in alpha
    for(int i = 0; i < n; i++) {
      double za = (Z[i] - 1) / (alpha * Z[i] + 1 - alpha);
      da += za;
      dda += za*za;
    }
    alpha1 = alpha + da/dda;
    alpha1 = (alpha1 < 1)?alpha1:1;
    if(alpha1 < 0) return NumericVector::create(NA_REAL, NA_REAL);
    // stopping criterion is |alpha - alpha1| < eps : we're not moving fast on the x axis
    // and da < 1 : we're not moving faster on the y axis than on the x axis
    // note that we tested the bounds 0 and 1 before entering the loop so we can't be at a bound with a large derivative
    if( std::abs(alpha - alpha1) < eps && da < 1 ) break;
    alpha = alpha1;
  }
  double hflod = 0;
  for(int i = 0; i < n; i++) hflod += std::log10(alpha1 * Z[i] + 1 - alpha1);
  return NumericVector::create(alpha1, hflod);
}
