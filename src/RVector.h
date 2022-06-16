// d'après la classe homonyme de RcppParallel

#ifndef _rvector__
#define _rvector__

#include <cstddef>

template <typename T>
class RVector {

   T * beg_data;
   T * end_data; 
public:   
   template <typename Source>
   inline explicit RVector(const Source & source) 
      : beg_data( const_cast<T*>(source.begin()) ), end_data(beg_data + source.length()) { }

   inline RVector(T * begin, T * end_) : beg_data(begin), end_data(end_) { }
   
   inline RVector(const RVector & other) : beg_data(other.beg_data), end_data(other.end_data) { }
   
   inline RVector & operator=(const RVector & rhs) {
      beg_data = rhs.beg_data;
      end_data = rhs.end_data;
      return *this;
   }

   inline RVector() : beg_data(nullptr), end_data(nullptr) {}
 
   inline T * begin() { return beg_data; }
   inline T * end() { return end_data; }
   
   inline const T * begin() const { return beg_data; }
   inline const T * end() const { return end_data; }
   
   inline std::size_t size() const { return end_data - beg_data; }
   inline std::size_t length() const { return end_data - beg_data; }
   
   inline T & operator[](std::size_t i) {
     return *(beg_data + i);
   }
   
   inline const T & operator[](std::size_t i) const {
     return *(beg_data + i);
   }
};

#endif 
