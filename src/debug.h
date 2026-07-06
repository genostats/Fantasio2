#include <iostream>

#ifndef SHOW
#define SHOW(x) std::cout << #x << " = " << (x) << "\n";
#endif

#ifndef SHOWVEC
#define SHOWVEC(x) { std::cout << #x << " = "; \
  for(auto & _a_ : x) std::cout << _a_ << " "; \
  std::cout << "\n"; }
#endif
