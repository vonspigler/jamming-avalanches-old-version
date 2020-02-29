#ifndef HEADERS_VECTORS
#define HEADERS_VECTORS

/*##################################################################################################
##                                                                                                ##
##   - class VECTORS, instantiated to DOUBLE and UNSIGNED INT (not needed anymore)                ##
##                                                                                                ##
##   It stores matrices into linear arrays; it can read/write.                                    ##
##                                                                                                ##
##################################################################################################*/

template <class T> class vectors {
private:
  T *data;

public:
  unsigned int NN, DD;
  vectors();
  vectors(unsigned int, unsigned int);
  vectors(const vectors &);
  ~vectors();

  vectors<T>& operator=(const vectors<T> &);

  inline T& operator()(unsigned int i, unsigned int d) {
//    if (i >= NN) std::cout << "ERROR, i (" << i << ") must be < N! (" << NN << ")" << std::endl;
//    if (d >= DD) std::cout << "ERROR, d (" << d << ") must be < DIM! (" << DD << ")" << std::endl;
    return data[i*DD + d];
  };
  inline const T& operator()(unsigned int i, unsigned int d) const {
//    if (i >= NN) std::cout << "ERROR, i (" << i << ") must be < N! (" << NN << ")" << std::endl;
//    if (d >= DD) std::cout << "ERROR, d (" << d << ") must be < DIM! (" << DD << ")" << std::endl;
    return data[i*DD + d];
  };
  inline const T* operator()(unsigned int i) const {
//    if (i >= NN) std::cout << "ERROR, i (" << i << ") must be < N! (" << NN << ")" << std::endl;
    return &data[i*DD];
  }; // READ-ONLY to avoid mess; call as const double* r = x(i);
};

extern template class vectors<double>;
extern template class vectors<unsigned int>;

#endif
