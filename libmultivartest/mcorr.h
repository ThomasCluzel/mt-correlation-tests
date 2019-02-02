//mcorr.h
#ifndef mcorr_H
#define mcorr_H

#include "./tntjama/tnt.h"
using namespace TNT;

class mcorr{
  //private class methods
  unsigned int mcorr_binomCoef(unsigned int N, 
			       unsigned int K);
  void mcorr_quickSort(TNT::Array2D<long double> matrix, 
		       int cols, int left, int right);

 public:
  //public class members
  int n; //number of rows of pseudorandom data
  int p; //number of columns of pseudorandom data
  TNT::Array2D <long double> mat;
  int numCorrs;
  double alpha;
  TNT::Array2D <long double> corrData;
  TNT::Array2D <long double> pVals;
  TNT::Array2D <long double> sortedPVals;
  TNT::Array2D <long double> Z;
  TNT::Array2D <long double> Ranks;

  //constructor & destructor
  mcorr(int N, int P, TNT::Array2D <long double> Mat,
	double Alpha);
  ~mcorr(){};

  //public class methods
  void mcorr_corr(TNT::Array2D <long double> matrix);
  int mcorr_BHY();
  void mcorr_fisherTrans(int type);
  void mcorr_getPVals();
  void mcorr_spearman();
  void mcorr_kendall();
  void mcorr_kendallNormal();
   //corrType -> 0 = Pearson, 1 = Spearman, 2 = Kendall
  void mcorr_pairCorr(int corrType);

};
#endif
