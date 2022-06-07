//[[Rcpp::export]]
bool checkOpenMP() {
#if defined(_OPENMP)
  return true;
#else
  return false;
#endif
}

