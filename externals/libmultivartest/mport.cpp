//mtimes.cpp
#include "./tntjama/jama_cholesky.h"
#include "./tntjama/jama_lu.h"
#include "./tntjama/jama_qr.h"
#include "./tntjama/tnt.h"
#include "mport.h"
#include <iomanip>
#include <cmath>
using std::cout;
using std::endl;
using std::ios;
using std::setw;

//Constructor
mport::mport(int N, int P, TNT::Array2D <long double> Mat, long double Alpha){
    cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl
	 << "                 Starting Multivariate Extension" << endl
	 << "                 Version: TestU01 1.2.3" << endl
	 << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl << endl << endl;

  if(N >= 1 && P >= 1){
    n = N;
    p = P;
    mat = Mat;
    alpha = Alpha;
    lagOrder = 20;  //Specify default value if none given
  }
  else
    std::cout << "Number of rows and columns must be positive integers." 
	      << std::endl;
};


//Private method definitions

//Creates identity matrix
TNT::Array2D<long double> mport::identity(int size){
  TNT::Array2D <long double> iMat(size, size, 0.0);
  for(int j = 0; j < iMat.dim1(); j++)
    iMat[j][j] = 1.0;
  return iMat;
}

//From http://wiki.cs.princeton.edu/index.php/TNT
//Compute transpose of a matrix
template<class T>
TNT::Array2D<T> mport::transpose(const TNT::Array2D<T> &M)
{
  TNT::Array2D<T> tran(M.dim2(), M.dim1() );
  for(int r = 0; r < M.dim1(); ++r)
    for(int c = 0; c < M.dim2(); ++c)
      tran[c][r] = M[r][c];
  return tran;
}

//Get diagonal values of matrix in 1D vector
template<class T>
TNT::Array1D<T> mport::diag(const TNT::Array2D<T> &M)
{
  TNT::Array1D<T> diagVec(M.dim1(), 0.0);
  for(int r = 0; r < M.dim1(); r++)
    diagVec[r] = M[r][r];
  return diagVec;
}


//Computes an outer product of two inputted vectors
template <class T>
TNT::Array2D<T> mport::mport_outerProd(const TNT::Array1D<T> &v, const TNT::Array1D<T> &v2)
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


//Added from http://wiki.cs.princeton.edu/index.php/TNT
//Invert matrix
TNT::Array2D <long double> mport::invert(const TNT::Array2D<long double> &M){
  assert(M.dim1() == M.dim2()); // square matrices only please

  // solve for inverse with LU decomposition
  JAMA::LU<long double> lu(M);

  // create identity matrix
  TNT::Array2D<long double> id(M.dim1(), M.dim2(), 0.0);
  for (int i = 0; i < M.dim1(); i++) 
    id[i][i] = 1;

  // solves A * A_inv = Identity
  return lu.solve(id);
}

//Mean center inputted matrix
void mport::mport_centerMat (){  
  TNT::Array2D <long double> meanMat;
  centMat = TNT::Array2D <long double> (mat.dim1(), mat.dim2(), 0.0);

  meanMat = TNT::Array2D <long double> (1, mat.dim2(), 0.0);
  //Means of columns of matrix
  for(int j = 0; j < mat.dim2(); j++){
    for(int i = 0; i < mat.dim1(); i++)
      meanMat[0][j] += mat[i][j];
    meanMat[0][j] /= mat.dim1();
  }

  TNT::Array2D <long double> repmat (mat.dim1(), mat.dim2(), 0.0);
  for(int i = 0; i < repmat.dim1(); i++)
    for(int j = 0; j < repmat.dim2(); j++)
      repmat[i][j] = meanMat[0][j];

  //Mean center
  centMat = mat - repmat;

}

//Used to construct Toeplitz Block matrix (matrix of matrices for GV test)
void mport::mport_blockFill(long double** pIn, TNT::Array2D <long double> pToFill,
		       int rowLoc, int colLoc, int blkSize){
  //put in the blkSize-by-blkSize array in the (rowLoc, colLoc) block
  if(colLoc > rowLoc) //top triangle
    for(int i = 0; i < blkSize; i++)
      for(int j = 0; j < blkSize; j++)
	pToFill[i + blkSize*rowLoc][j + blkSize*colLoc] = pIn[i][j];
  else //otherwise put in the transpose (bottom triangle)
    for(int i = 0; i < blkSize; i++)
      for(int j = 0; j < blkSize; j++)
	pToFill[i + blkSize*rowLoc][j + blkSize*colLoc] = pIn[j][i];
}

//Hosking and Li-McLeod tests
void mport::mport_portmanteauTests(int lagOrder){
  TNT::Array2D <long double> c0(p, p, 0.0);
  TNT::Array2D <long double> c0inv(c0.dim1(), c0.dim2(), 0.0);
  TNT::Array1D <long double> d;
  TNT::Array2D <long double> dd;
  TNT::Array2D <long double> L;
  TNT::Array2D <long double> matcTemp1;
  TNT::Array2D <long double> matcTemp2;
  TNT::Array2D <long double> tempMult;
  TNT::Array3D <long double> cl(lagOrder, p, p, 0.0);
  long double lmp = 0.0; //Li-McLeod Portmanteau observed test statistic
  long double impCalc = 0.0; //Vectorized/Kronecker product used in Hosking and LM
  long double tempDiv = 0.0;
  long double hoskTemp = 0.0;
  TNT::Array2D <long double> vec_cl(p * p, 1, 0.0);
  TNT::Array2D <long double> innerMult;
  double pValue_lmp = 0.0;
  long double hosk = 0.0; //Hosking Portmanteau observed test statistic
  double pValue_hosk = 0.0;

  //Compute lag zero correlation matrix
  c0 = matmult(transpose(centMat), centMat);
  d = diag(c0);
  dd = mport_outerProd(d, d);
  
  for(int i = 0; i < dd.dim1(); i++)
    for(int j = 0; j < dd.dim2(); j++)
      dd[i][j] = sqrt(dd[i][j]);
  
  for(int i = 0; i < dd.dim1(); i++)
    for(int j = 0; j < dd.dim2(); j++)
      c0[i][j] /= dd[i][j];
  
  c0inv = invert(c0);

  //Compute lag el (el = 0, ... lagOrder - 1) correlation matrices and store in 3D array
  for(int el = 0; el < cl.dim1(); el++){
    matcTemp1 = TNT::Array2D <long double> (n - el - 1, centMat.dim2(), 0.0);
    matcTemp2 = TNT::Array2D <long double> (n - el - 1, centMat.dim2(), 0.0);
    for(int r = 0; r < matcTemp1.dim1(); r++)
      for(int c = 0; c < matcTemp1.dim2(); c++)
	matcTemp1[r][c] = centMat[r][c];

    for(int r = 0; r < matcTemp2.dim1(); r++)
      for(int c = 0; c < matcTemp2.dim2(); c++)
	matcTemp2[r][c] = centMat[r + el + 1][c];

    tempMult = matmult(transpose(matcTemp1), matcTemp2);
    for(int r = 0; r < cl.dim2(); r++)
      for(int c = 0; c < cl.dim3(); c++)
	cl[el][r][c] = tempMult[r][c];

    for(int i = 0; i < dd.dim1(); i++)
      for(int j = 0; j < dd.dim2(); j++)
	cl[el][i][j] /= dd[i][j];    
  }

  //Compute covariance matrix in Hosking & Li-McLeod statistics
  
  //Compute Kronecker product of c0inv with itself
  TNT::Array2D <long double> rr(c0inv.dim1() * c0inv.dim1(), 
				c0inv.dim2() * c0inv.dim2(), 0.0);
  for(int r = 0; r < c0inv.dim1(); r++)
    for(int c = 0; c < c0inv.dim2(); c++)
      for(int i = 0; i < c0inv.dim1(); i++)
	for(int j = 0; j < c0inv.dim2(); j++)
	  rr[r * c0inv.dim1() + i][c * c0inv.dim2() + j] = c0inv[r][c] * c0inv[i][j];

  //Vectorize cl matrices
  for(int el = 0; el < cl.dim1(); el++){
    for(int i = 0; i < cl.dim2(); i++)
      for(int j = 0; j < cl.dim3(); j++)
	vec_cl[i + j * cl.dim2()][0] = cl[el][i][j];
    
    innerMult = matmult(matmult(transpose(vec_cl), rr), vec_cl);
    impCalc += innerMult[0][0];
    tempDiv = innerMult[0][0] / (long double) (n - el - 1);
    hoskTemp += tempDiv;
  }

  //Li-McLeod test results
  lmp = n * impCalc + p * p * lagOrder * (lagOrder + 1) / (long double) (2.0 * n);
  if(std::isnan(lmp)){
    cout << "Li-McLeod Test Statistic is infinite.  Program exiting." << endl;
    return;
  }
  df = p * p * lagOrder; //modelOrder = 0
  pValue_lmp = 1 - gamma_inc( (double) df / 2.0,  (double) lmp / 2.0);

  cout << "-----------------------------------------------" << endl
       << "Portmanteau Test for White Noise  (Li-McLeod):" << endl
       << "Alpha = " << alpha << endl << endl
       << "-----------------------------------------------" << endl
       << "   n =  " << n << ",  p = " << p << ",  Lag order = " << lagOrder
       << endl << endl << endl;

 cout << "-----------------------------------------------" << endl
      << "Li-McLeod Test Statistic            :  " << lmp << endl
      << "Degrees of Freedom                  :  " << df << endl
      << "p-value of test                     :  " << pValue_lmp << endl << endl << endl;

  /*
  if(pValue_lmp <= alpha)
    cout << "For the Li-McLeod test, since the p-value=" << pValue_lmp << " <= alpha="
	 << alpha << ", we reject the null hypothesis that the inputted pseudorandom "
	 << "vectors are white noise.  We suspect the vectors are dependent."
	 << endl << endl;
  else
    cout << "For the Li-McLeod test, there is not sufficient evidence to conclude"
	 << " that the residual vectors are correlated."
	 << endl << endl;
  */

  
  //Hosking test results
  hosk = n * n * hoskTemp;
  if(std::isnan(hosk)){
    cout << "Hosking Test Statistic is infinite.  Program exiting." << endl;
    return;
  }
  pValue_hosk = 1 - gamma_inc ( (double) df / 2.0, (double) hosk / 2.0);
  
  cout << "-----------------------------------------------" << endl
       << "Portmanteau Test for White Noise  (Hosking):" << endl
       << "Alpha = " << alpha << endl << endl
       << "-----------------------------------------------" << endl
       << "   n =  " << n << ",  p = " << p << ",  Lag order = " << lagOrder
       << endl << endl << endl;

 cout << "-----------------------------------------------" << endl
      << "Hosking Test Statistic              :  " << hosk << endl
      << "Degrees of Freedom                  :  " << df << endl
      << "p-value of test                     :  " << pValue_hosk << endl << endl << endl;


  /*
  if(pValue_hosk <= alpha)
    cout << "For the Hosking test, the p-value=" << pValue_hosk << " <= alpha="
	 << alpha << ", we reject the null hypothesis that the inputted pseudorandom"
	 << "vectors are white noise.  We suspect the vectors are dependent."
	 << endl << endl;
  else
    cout << "For the Hosking test, there is not sufficient evidence to conclude"
	 << " that the inputted pseudorandom vectors are correlated."
	 << endl << endl;
  */
}

//Mahdi and McLeod's Generalized Variance test
void mport::mport_mahdiMcLeod(int lagOrder){
  TNT::Array2D <long double> cov0(p, p, 0.0);
  TNT::Array2D <long double> cov0inv(cov0.dim1(), cov0.dim2(), 0.0);
  TNT::Array2D <long double> matcTemp1;
  TNT::Array2D <long double> matcTemp2;
  TNT::Array2D <long double> tempMult;
  TNT::Array3D <long double> covl(lagOrder, p, p, 0.0);
  TNT::Array3D <long double> Rl(lagOrder, p, p, 0.0);
  TNT::Array2D <long double> mahdiMat(Rl.dim2() * (lagOrder + 1), Rl.dim3() * (lagOrder + 1), 0.0);
  long double gv = 0.0;
  long double df_gv = 0.0;
  double pValue_gv = 0.0;

  //Compute lag zero covariance matrix
  cov0 = matmult(transpose(centMat), centMat);

  for(int i = 0; i < cov0.dim1(); i++)
    for(int j = 0; j < cov0.dim2(); j++)
      cov0[i][j] /= (long double) n;
  
  // LL' = cov0^{-1}
  cov0inv = invert(cov0);
  JAMA::Cholesky<long double> chol(cov0inv);
  TNT::Array2D <long double> L = chol.getL();

  //Compute lag el covariance matrix
  for(int el = 0; el < covl.dim1(); el++){
    matcTemp1 = TNT::Array2D <long double> (n - el - 1, centMat.dim2(), 0.0);
    matcTemp2 =  TNT::Array2D <long double> (n - el - 1, centMat.dim2(), 0.0);
    for(int r = 0; r < matcTemp1.dim1(); r++)
      for(int c = 0; c < matcTemp1.dim2(); c++)
	matcTemp1[r][c] = centMat[r][c];

    for(int r = 0; r < matcTemp2.dim1(); r++)
      for(int c = 0; c < matcTemp2.dim2(); c++)
	matcTemp2[r][c] = centMat[r + el + 1][c];

    tempMult = matmult(transpose(matcTemp1), matcTemp2);
    for(int r = 0; r < covl.dim2(); r++)
      for(int c = 0; c < covl.dim3(); c++)
	covl[el][r][c] = tempMult[r][c];

    
    for(int i = 0; i < covl.dim2(); i++)
      for(int j = 0; j < covl.dim3(); j++)
	covl[el][i][j] /= (long double) n;
  }

  //Rl = L' * covl * L (matrix multiplication)
  for(int el = 0; el < lagOrder; el++){
    TNT::Array2D <long double> tempL(p, p, 0.0);
    TNT::Array2D <long double> tempL2(p, p, 0.0);
    for(int i = 0; i < Rl.dim2(); i++)
      for(int j = 0; j < Rl.dim3(); j++)
	tempL[i][j] = covl[el][i][j];
    tempL2 = transpose(matmult(matmult(transpose(L), tempL), L));
    for(int i = 0; i < Rl.dim2(); i++)
      for(int j = 0; j < Rl.dim3(); j++)
	Rl[el][i][j] = tempL2[i][j];
  }

  //Mahdi and McLeod Generalized Variance test

  //Initialize diagonal to 1
  for(int i = 0; i < mahdiMat.dim1(); i++)
    for(int j = 0; j < mahdiMat.dim2(); j++){
      if(i == j)
	mahdiMat[i][j] = 1.0;
      else 
	mahdiMat[i][j] = 0.0;
    }

  //Fill the upper triangle of the matrix of matrices
  for(int i = 0; i < lagOrder + 1; i++)
    for(int j = i + 1; j < lagOrder + 1; j++)
      mport_blockFill(Rl[j - i - 1], mahdiMat, i, j, Rl.dim2());

  //now fill the lower part
  for(int i = 1; i < lagOrder + 1; i++)
    for(int j = 0; j < i; j++)
      mport_blockFill(Rl[i - j - 1], mahdiMat, i, j, Rl.dim2());

  JAMA::LU<long double> lu(mahdiMat);

  //Generalized Variance test results
  gv = (-3.0 * n) / (2.0 * lagOrder + 1.0) * log(lu.det());
  if(std::isnan(gv)){
    cout << "GV Test Statistic is infinite.  Program exiting." << endl;
    return;
  }
  df_gv = p * p * (1.5 * lagOrder * (lagOrder + 1) / (2.0 * lagOrder + 1.0) ); //modelOrder = 0
  pValue_gv = 1 - gamma_inc ( (double) df_gv / 2.0, (double) gv / 2.0);

  cout << "-----------------------------------------------" << endl
       << "Portmanteau Test for White Noise  " << endl
       << "   (Mahdi-McLeod):" << endl
       << "Alpha = " << alpha << endl << endl
       << "-----------------------------------------------" << endl
       << "   n =  " << n << ",  p = " << p << ",  Lag order = " << lagOrder
       << endl << endl << endl;

 cout << "-----------------------------------------------" << endl
      << "Determinant of Toeplitz Block matrix:  " << lu.det() << endl
      << "Mahdi-McLeod Test Statistic         :  " << gv << endl
      << "Degrees of Freedom                  :  " << df_gv << endl
      << "p-value of test                     :  " << pValue_gv << endl << endl << endl;


  /*
  if(pValue_gv <= alpha)
    cout << "For the Generalized Variance Test, since the p-value=" << pValue_gv
	 << " <= alpha=" << alpha 
	 << ", we reject the null hypothesis that the inputted pseudorandom"
	 << "vectors are white noise.  We suspect the vectors are dependent."
	 << endl << endl;
  else
    cout << "For the GV test, there is not sufficient evidence to conclude"
	 << " that the inputted pseudorandom vectors are correlated."
	 << endl << endl;
  */
}

