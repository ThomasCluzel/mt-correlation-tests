//mport.h
#ifndef mport_H
#define mport_H

#include "./prob/prob.hpp"

class mport{
  //private class methods
  void mport_blockFill(long double** pIn, TNT::Array2D <long double> pToFill,
		 int rowLoc, int colLoc, int blkSize);
  TNT::Array2D <long double> identity(int size);
  template<class T> TNT::Array2D<T> transpose(const TNT::Array2D<T> &M);
  template<class T> TNT::Array1D<T> diag(const TNT::Array2D<T> &M);
  template<class T> TNT::Array2D<T> mport_outerProd(const TNT::Array1D<T> &v, const TNT::Array1D<T> &v2);
  TNT::Array2D <long double> invert(const TNT::Array2D<long double> &M);

 public:
  //public class members
  int n; //number of rows of pseudorandom data
  int p; //number of columns of pseudorandom data
  TNT::Array2D <long double> mat; //inputted matrix of pseudorandom data
  TNT::Array2D <long double> centMat; //mean centered matrix to test for white noise
  long double alpha;  //to compare p-values against
  int lagOrder; //upper limit of summation in test statistics
  int df; //degrees of freedom for Hosking & Li-McLeod tests
  int mahdiDf; //degrees of freedom for Mahdi-McLeod Generalized Variance tests
  
  //constructor & destructor
  mport(int N, int P, TNT::Array2D <long double> Mat, long double Alpha);
  ~mport(){};

  //public class methods
  void mport_centerMat(); //mean center inputted matrix values
  void mport_portmanteauTests(int lagOrder); //Li-McLeod and Hosking Tests
  void mport_mahdiMcLeod(int lagOrder); //Generalized Variance test

};
#endif
