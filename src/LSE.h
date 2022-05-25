
// LSE(a, b) = log( e^a + e^b )
#define LSE(a,b) (((a)>(b))?(a)+log1p(exp((b)-(a))):(b)+log1p(exp((a)-(b))))

/* ceci m'envoie des warnings à la pelle en -pedantic
  #define LSE(a,b) ({ \
    __typeof__(a) __a=(a); \
    __typeof__(b) __b=(b); \
    (__a>__b)?__a+log1p(exp(__b-__a)):__b+log1p(exp(__a-__b)); \
  })
*/

/*
template<typename T>
inline T LSE(T a, T b) {
  return (a>b)?a+log1p(exp(b-a)):b+log1p(exp(a-b));
}
*/
