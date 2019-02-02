//mmult.h
#ifndef mmult_H
#define mmult_H

#include "./prob/prob.hpp"

class mmult{
  //private class methods
  long double trace (TNT::Array2D <long double> A);
  TNT::Array2D <long double> multConst (TNT::Array2D <long double> A, long double b);
  template<class T> TNT::Array1D<T> copyRowToVec(const TNT::Array2D<T> &M, int rowNum);
  template<class T> TNT::Array2D<T> outerProd(const TNT::Array1D<T> &v, 
					      const TNT::Array1D<T> &v2);

 public:
  //public class members
  int n; //number of rows of pseudorandom data
  int p; //number of columns of pseudorandom data
  TNT::Array2D <long double> mat; //inputted matrix
  double alpha; //significance level
  TNT::Array2D <long double> normMat;  //normalized input
  TNT::Array2D <long double> A; //related to covariance matrix
  TNT::Array2D <long double> C; //n^{-1} A
  long double obsTS; //observed likelihood ratio test statistic
  double pValue; //corresponding p-value

  //constructor & destructor
  mmult(int N, int P, TNT::Array2D <long double> Mat, double Alpha);
  ~mmult(){};

  //public class methods
  void mmult_LRT();

};
#endif
