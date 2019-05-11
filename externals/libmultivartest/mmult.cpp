//mmult.cpp
#include "./tntjama/tnt.h"
#include "./tntjama/jama_lu.h"
using namespace TNT;
#include "mmult.h"
#include <cmath>
using std::cout;
using std::endl;
using std::ios;
using std::swap;

//Constructor
mmult::mmult(int N, int P, TNT::Array2D <long double> Mat, double Alpha ){
    cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl
	 << "                 Starting Multivariate Extension" << endl
	 << "                 Version: TestU01 1.2.3" << endl
	 << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl << endl << endl;

  if(N >= 1 && P >= 1){
    n = N;
    p = P;
    mat = Mat;
    alpha = Alpha;
    normMat = TNT::Array2D <long double> (N, P, 0.0);
    A = TNT::Array2D <long double> (P, P, 0.0);
    C = TNT::Array2D <long double> (A.dim1(), A.dim2(), 0.0);
    obsTS = 0.0;
    pValue = 0.0;
  }
  else
    std::cout << "Number of rows and columns must be positive integers." 
	      << std::endl;
};

//Private method definitions
//Determine trace of inputted square matrix
long double mmult::trace(TNT::Array2D <long double> A){
  long double matTrace = 0.0;
  for(int i = 0; i < A.dim1(); i++)
    matTrace += A[i][i];
  return matTrace;
}

//Multiply matrix by a constant
TNT::Array2D <long double> mmult::multConst (TNT::Array2D <long double> A, long double b){
  TNT::Array2D <long double> newMat (A.dim1(), A.dim2(), 0.0);
  for(int i = 0; i < newMat.dim1(); i++)
    for(int j = 0; j < newMat.dim2(); j++)
      newMat[i][j] = b * A[i][j];

  return newMat;
}

//Copy matrix row into Array1D
template<class T>
TNT::Array1D<T> mmult::copyRowToVec(const TNT::Array2D<T> &M, int rowNum){
  TNT::Array1D<T> vec(M.dim2());
  for(int c = 0; c < M.dim2(); ++c)
    vec[c] = M[rowNum][c];
  return vec;
}

//Computes an outer product of two inputted vectors
template <class T>
TNT::Array2D<T> mmult::outerProd(const TNT::Array1D<T> &v, const TNT::Array1D<T> &v2)
{    
  //declare variable to store matrix
  TNT::Array2D<T> outerMat(v.dim(), v2.dim(), 0.0);
    
  //multiply components in vector
  for(int i = 0; i < v.dim(); i++)
    for(int j = 0; j < v2.dim(); j++)
      outerMat[i][j] = v[i] * v2[j];
    
  //return answer
  return outerMat;
}

//Public method definitions
void mmult::mmult_LRT(){

  //Convert uniform (0, 1) deviate matrix to normal(0, 1) deviates
  for(int i = 0; i < mat.dim1(); i++)
    for(int j = 0; j < mat.dim2(); j++)
      normMat[i][j] = normal_01_cdf_inv (mat[i][j]); //from prob.cpp

  //Compute A matrix 
  for(int i = 0; i < n; i++){
    TNT::Array1D <long double> temp1(normMat.dim2(), 0.0);
    temp1 = copyRowToVec(normMat, i);
    A += outerProd(temp1, temp1);
  }

  //Muirhead calculations using C = n^{-1} A
  C = multConst(A,  (long double) (1.0 /  (long double) n));

  JAMA::LU<long double> luC(C);
  
  obsTS = n * (trace(C) - log(luC.det()) - p);

  cout << "-----------------------------------------------" << endl
       << "Likelihood Ratio Test for " << endl
       << "     Pairwise Correlation Matrix = Identity:" << endl
       << "Alpha = " << alpha << endl << endl
       << "-----------------------------------------------" << endl
       << "   n =  " << n << ",  p = " << p << endl << endl << endl;

  cout << "-----------------------------------------------" << endl
       << "LR Test Statistic                   :  " << obsTS << endl
       << "p-value of test                     :  ";

  if(!std::isnan(obsTS))
    pValue = 1 - chi_square_cdf(obsTS, p * (p + 1.0) / 2.0) ;
  cout << pValue << endl << endl;


}
