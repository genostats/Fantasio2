#ifndef _RMATRIX_GASTON2_
#define _RMATRIX_GASTON2_

#include <cstddef>

template <typename T>
class RMatrix { 
    
    private : 
    T * beg_data;
    T * end_data; // used for bound checking
    size_t nrow_; // for double indexes access
    size_t ncol_; // ncol is deductible theoretically

    public :
    using value_type = T;
    
    inline T * begin() { return beg_data; }
    inline T * end() { return end_data; }

    inline const T * begin() const { return beg_data; }
    inline const T * end() const { return end_data; }


    inline std::size_t nrow() const { return nrow_; }
    inline std::size_t ncol() const { return ncol_; }


    inline std::size_t size() const { return nrow_ * ncol_; }
    inline std::size_t length() const { return nrow_ * ncol_; } // or end_data - beg_data ?


    // ----------------- constructors --------------------------

    inline RMatrix() : beg_data(nullptr), end_data(nullptr), nrow_(0), ncol_(0) {}
    inline RMatrix(T * begin, T * end_, size_t nrow, size_t ncol) : beg_data(begin), end_data(end_),  nrow_(nrow), ncol_(ncol) { }
    // constructeur depuis n'importe quel autre type de datastruct qui supporte .begin(), .length(), .nrow(), .ncol()
    template <typename Source>
    inline explicit RMatrix(const Source & source) 
      : beg_data( const_cast<T*>(source.begin()) ), end_data(beg_data + source.length()), nrow_(source.nrow()), ncol_(source.ncol()) { }
     
    inline RMatrix(const RMatrix & other) : beg_data(other.beg_data), end_data(other.end_data), nrow_(other.nrow()), ncol_(other.ncol()) { }
   

    // ----------------- operator [] --------------------------
    // Operator [] gives back the data at index, UNSAFE.
    inline T & operator[](size_t ind) {
    return beg_data[ind]; // *(beg_data + ind);
    }

    inline const T & operator[](size_t ind) const
    {
    return beg_data[ind]; //*(beg_data + ind);
    }

    // ----------------- operator () --------------------------
    // data at row i, col j
    inline T & operator()(size_t i, size_t j)
    {
    return beg_data[(j * nrow_) + i];
    }

    inline const T & operator()(size_t i, size_t j) const
    {
    return beg_data[(j * nrow_) + i];
    }

    // ----------------- operator = --------------------------

    inline RMatrix & operator=(const RMatrix & rhs) {
      beg_data = rhs.beg_data;
      end_data = rhs.end_data;
      nrow_ = rhs.nrow_;
      ncol_ = rhs.ncol_;
      return *this;
   }
};

#endif //_RMATRIX_GASTON2_