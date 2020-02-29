#include                    "vectors.h"
#include                    <cstdlib>
#include                    <cstring>
#include                    <iostream>

template <class T> vectors<T>::vectors(const vectors<T> &copy) {
  data = new T [copy.NN*copy.DD];
  std::memcpy(data, copy.data, copy.NN*copy.DD*sizeof(T));
  NN = copy.NN;
  DD = copy.DD;
}

template <class T> vectors<T>::vectors() : NN(-1), DD(-1) { data = new T [1]; }

template <class T> vectors<T>::vectors(unsigned int N, unsigned int D) : NN(N), DD(D) {
  data = new T [N*D];
}

template <class T> vectors<T>::~vectors() { delete [] data; }

template <class T> vectors<T>& vectors<T>::operator= (const vectors<T> &copy) {
  if (this != &copy) {
    delete [] data;
    data = new T [copy.NN*copy.DD];
    std::memcpy(data, copy.data, copy.NN*copy.DD*sizeof(T));
    NN = copy.NN;
    DD = copy.DD;
  }
  return *this;
}
/*template <class T> T& vectors<T>::operator()(unsigned int i, unsigned int d) {
  if (i >= NN) std::cout << "ERROR, i (" << i << ") must be < N! (" << NN << ")" << std::endl;
  if (d >= DD) std::cout << "ERROR, d (" << d << ") must be < DIM! (" << DD << ")" << std::endl;
  return data[i*DD + d];
}

template <class T> const T& vectors<T>::operator()(unsigned int i, unsigned int d) const {
  if (i >= NN) std::cout << "ERROR, i (" << i << ") must be < N! (" << NN << ")" << std::endl;
  if (d >= DD) std::cout << "ERROR, d (" << d << ") must be < DIM! (" << DD << ")" << std::endl;
  return data[i*DD + d];
}

template <class T> const T* vectors<T>::operator()(unsigned int i) const {                        // READ-ONLY to avoid mess
  if (i >= NN) std::cout << "ERROR, i (" << i << ") must be < N! (" << NN << ")" << std::endl;
  return &data[i*DD];
}*/

template vectors<double>::vectors();
template vectors<double>::vectors(unsigned int, unsigned int);
template vectors<double>::vectors(const vectors<double> &);
template vectors<double>::~vectors();
template vectors<double>& vectors<double>::operator= (const vectors<double> &);
template double& vectors<double>::operator()(unsigned int, unsigned int);
template const double& vectors<double>::operator()(unsigned int, unsigned int) const;
template const double* vectors<double>::operator()(unsigned int) const;

template vectors<unsigned int>::vectors();
template vectors<unsigned int>::vectors(unsigned int, unsigned int);
template vectors<unsigned int>::vectors(const vectors<unsigned> &);
template vectors<unsigned int>::~vectors();
template vectors<unsigned int>& vectors<unsigned int>::operator=(const vectors<unsigned int> &);
template unsigned int& vectors<unsigned int>::operator()(unsigned int, unsigned int);
template const unsigned int& vectors<unsigned int>::operator()(unsigned int, unsigned int) const;
template const unsigned int* vectors<unsigned int>::operator()(unsigned int) const;
