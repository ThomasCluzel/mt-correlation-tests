//mcorr.cpp
#include "./tntjama/tnt.h"
using namespace TNT;
#include "kendall2.c"
#include "mcorr.h"
#include <cmath>
using std::cout;
using std::endl;
using std::ios;
using std::swap;

//Constructor
mcorr::mcorr(int N, int P, TNT::Array2D <long double> Mat, double Alpha ){
    cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl
	 << "                 Starting Multivariate Extension" << endl
	 << "                 Version: TestU01 1.2.3" << endl
	 << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl << endl << endl;

  if(N >= 1 && P >= 1){
    n = N;
    p = P;
    mat = Mat;
    numCorrs = mcorr_binomCoef(p, 2);
    alpha = Alpha;
    corrData = TNT::Array2D <long double>(numCorrs, 3);
    pVals = TNT::Array2D <long double>(numCorrs, 3);
    Z = TNT::Array2D <long double>(numCorrs, 3);
    Ranks = TNT::Array2D <long double>(n, p);
  }
  else
    std::cout << "Number of rows and columns must be positive integers." 
	      << std::endl;
};

//Private method definitions
unsigned int mcorr::mcorr_binomCoef(unsigned int N, unsigned int K){
  if(K == 0 || K == N) 
    return 1;
  else 
    return (N*mcorr_binomCoef(N - 1, K - 1))/K;
}

void mcorr::mcorr_quickSort(TNT::Array2D<long double> matrix, int cols, int left, int right) {
      int i = left, j = right;
      long double pivot = matrix[(left + right) / 2][0];
 
      // Partitioning
      while (i <= j) {
	while (matrix[i][0] < pivot)
	  i++;
	while (matrix[j][0] > pivot)
	  j--;
	if (i <= j) {
	  std::swap(matrix[i][0], matrix[j][0]);
	  for(int k = 1; k < cols; k++){
	    std::swap(matrix[i][k], matrix[j][k]);
	  }
	  i++;
	  j--;
	}
      }
 
      // Recursive calls
      if (left < j)
	mcorr_quickSort(matrix, cols, left, j);
      if (i < right)
	mcorr_quickSort(matrix, cols, i, right);
}

//Public method definitions
void mcorr::mcorr_corr(TNT::Array2D <long double> matrix){
    int count = 0;
    TNT::Array1D <long double> mean(p, 0.0);
    TNT::Array2D <long double> Cov(p, p, 0.0);

     // Two pass algorithm
    //Calculate the means of each column
    for(int j = 0; j < p; j++){
      for (int i = 0; i < n; i++)
    	mean[j] += mat[i][j];
      mean[j] /= n;
    }
    //Compute covariance matrix
    for(int k = 0; k < p; k++){
      for(int j = 0; j < p; j++){
  	for (int i = 0; i < n; i++)
  	  Cov[k][j] += (mcorr::mat[i][j] - mean[j])*(mcorr::mat[i][k] - mean[k]);
  	Cov[k][j] /= (n - 1);
      }
    }
    //  cout << "Covariance:  " << Cov;

    //Compute the Pearson correlation coefficient matrix 
      //(only bottom triangular to save space and time)
    for(int k = 0; k < mcorr::p; k++)
      for(int j = 0; j < k; j++){
  	if(Cov[k][k]*Cov[j][j] == 0.0 && Cov[k][j] >= 0)
	  mcorr::corrData[count][0] = 1.0000;
  	else if(Cov[k][k]*Cov[j][j] == 0.0 && Cov[k][j] < 0)
	  mcorr::corrData[count][0] = -1.0000;
  	else
	  mcorr::corrData[count][0] = Cov[k][j]/sqrt(Cov[k][k]*Cov[j][j]); 
    
	mcorr::corrData[count][1] = k;
	mcorr::corrData[count][2] = j;
	count++;
      }
}

void mcorr::mcorr_fisherTrans(int type){
  for(int i = 0; i < mcorr::numCorrs; i++){
    mcorr::Z[i][0] = sqrt( (mcorr::n - 3)/1.06 ) * 0.5 * log( (1.0 + mcorr::corrData[i][0]) / (1.0 - mcorr::corrData[i][0]) );
    if (type == 1)
      mcorr::Z[i][0] *= sqrt(1.06);
    mcorr::Z[i][1] = mcorr::corrData[i][1];
    mcorr::Z[i][2] = mcorr::corrData[i][2];
  }
}

void mcorr::mcorr_getPVals(){
  for(int i = 0; i < mcorr::numCorrs; i++){
    //Two tailed P-value from Z test
    mcorr::pVals[i][0] = erfc(fabs(mcorr::Z[i][0]) / sqrt(2) );
    mcorr::pVals[i][1] = mcorr::Z[i][1];
    mcorr::pVals[i][2] = mcorr::Z[i][2];
  }
}

int mcorr::mcorr_BHY(){
  double q;
  int k;
  sortedPVals = TNT::Array2D <long double> (mcorr::numCorrs, 4);
  double tempSum;
  
  //Arrange P-values in ascending order
  mcorr_quickSort(mcorr::pVals, 3, 0, mcorr::numCorrs - 1);
  /*  Simple sort for check
    for(int i = 1; i < m; i++)
    for(int j = 1; j < m; j++)
      if(pVals[j][0] < pVals[j-1][0]){
	std::swap(pVals[j][0], pVals[j-1][0]);
	std::swap(pVals[j][1], pVals[j-1][1]);
	std::swap(pVals[j][2], pVals[j-1][2]);
      }
  */

  //Prepare rankings
  for(int i = 0; i < mcorr::numCorrs; i++){
    sortedPVals[i][0] = i;
    for(int j = 0; j < 3; j++)
      sortedPVals[i][j+1] = mcorr::pVals[i][j];
  }

  //Define q
  for(double j = 1.0; j <= mcorr::numCorrs; j++)
    tempSum += 1/j;
  q = alpha/tempSum;

  //Compute k
  k = 0;
  //cout << "Sorted P-Values Comparison" << endl;
  for(double i = 1.0; i <= mcorr::numCorrs; i++){
    /*  cout << sortedPVals[(int)(i-1)][1] << "--" 
	 << q << "*" << (int)i << "/" << numCorrs << "=" 
	 << q*(double)(i/mcorr::numCorrs) << endl;
    */
    if(sortedPVals[i-1][1] <= q*(long double)(i/mcorr::numCorrs))
      k++;
  }
  //cout << endl;

  return k;
  /*
  if(k == 0){
    cout << "Reject none of the null hypotheses" << endl;
    return;
  }
  
  cout << "Reject the null hypotheses of nonzero correlation corresponding to" << endl;
  for(int i = 0; i < k; i++){
    cout << "P-value_(" << i + 1 << "):  Vector " << (int)(sortedPVals[i][3] + 1) 
	 << " and Vector " << (int)(sortedPVals[i][2] + 1) << endl;
  }
  */
}


void mcorr::mcorr_spearman(){ 
  TNT::Array1D<long double> rankings(n, 0.0);
  TNT::Array2D<long double> temp(n, 3, 0.0);
  int j, ji, jt;
  double rank;

  for(int k = 0; k < p; k++){
    for(int i = 0; i < n; i++){
      temp[i][0] = mcorr::mat[i][k];
      temp[i][1] = i; //Original position (to sort by later)
    }
    j = 1;

    //Sort by first column and copy sorted list to rankings
    mcorr_quickSort(temp, 2, 0, n - 1);
    for(int i = 0; i < n; i++)
      rankings[i] = temp[i][0];


    //Rank the sorted vector (including midranks for ties)
    while (j < n) {
      if (rankings[j] != rankings[j - 1]) { //Not a tie.
	rankings[j-1] = j;
	++j;
      } 
      else { //A tie:
	for (jt = j + 1; jt <= n && rankings[jt - 1] == rankings[j - 1]; jt++);
	rank = 0.5 * (j + jt - 1); //Mean rank of the tie
	for (ji = j; ji <= (jt - 1); ji++) //Enter mean rank into all tied entries
	  rankings[ji - 1] = rank;
	j = jt;
      }
    }
    //If the last element was not tied, this is its rank
    if (j == n) 
      rankings[n - 1] = n; 

    //Swap first column and second column and insert rankings as third column of temp
    for(int i = 0; i < n; i++){
      temp[i][2] = rankings[i];
      swap(temp[i][0], temp[i][1]);
    }

    //Sort by original position
    mcorr_quickSort(temp, 3, 0, n - 1);

    //Place each rank column into the rank matrix that will be passed
    //into the Pearson correlation function
    for(int i = 0; i < n; i++)
      mcorr::Ranks[i][k] = temp[i][2];
  }

  //Compute Spearman correlation matrix (only bottom triangular)
  mcorr::mcorr_corr(Ranks);//corr(Ranks, corrData, n, p); 
   
}


void mcorr::mcorr_kendall(){
  //Since kendallNlogN requires double*'s
  double* arr1 = new double[n];
  double* arr2 = new double[n];

  int m = 0;
  TNT::Array2D <long double> newtemp(n, 2, 0.0);
   for(int k = 0; k < p; k++)
      for(int j = 0; j < k; j++){

	for(int i = 0; i < n; i++){
	  newtemp[i][0] = mcorr::mat[i][k];  
	  newtemp[i][1] = mcorr::mat[i][j];
	}

	//Sort in lockstep by column 1
	mcorr_quickSort(newtemp, 2, 0, n - 1);

	for(int i = 0; i < n; i++){
	  arr1[i] = newtemp[i][0];  
	  arr2[i] = newtemp[i][1];
	}
    
	mcorr::corrData[m][0] = kendallNlogN(arr1, arr2, n, 1);
	mcorr::corrData[m][1] = k;
	mcorr::corrData[m][2] = j;
	m++;
      }
  delete[] arr1;
  delete[] arr2;
}

void mcorr::mcorr_kendallNormal(){
  for(int i = 0; i < mcorr::numCorrs; i++){
    mcorr::Z[i][0] = mcorr::corrData[i][0] / sqrt( (2.0* (2.0*n + 5.0)) / (9.0 * n * (n - 1) ) );
    mcorr::Z[i][1] = mcorr::corrData[i][1];
    mcorr::Z[i][2] = mcorr::corrData[i][2];
  }
}

void mcorr::mcorr_pairCorr(int corrType){
  int k = 0;

  cout << "-----------------------------------------------" << endl
       << "PairCorr test ";
  if (corrType == 0)
    cout << "for Pearson Correlations:" << endl;
  else if (corrType == 1)
    cout << "for Spearman Correlations:" << endl;
  else
    cout << "for Kendall Correlations:" << endl;
  cout << "-----------------------------------------------" << endl
       << "   n =  " << n << ",  p = " << p << endl << endl << endl;
  if (corrType == 0)
    mcorr::mcorr_corr(mat);
  else if (corrType == 1)
    mcorr::mcorr_spearman();
  else
    mcorr::mcorr_kendall();
  if (corrType == 0){
    //Transform Pearson correlations into normals using the Fisher r to z Transform
    mcorr::mcorr_fisherTrans(1);
  }
  else if (corrType == 1){
    //Transform Spearman correlations into normals using the Fisher r to z Transform
    mcorr::mcorr_fisherTrans(0);
  }
  else{
    //Transform Kendall correlations into normals
    mcorr::mcorr_kendallNormal();
  }
  /*  cout << "-----------------------------------------------" << endl;
  cout << "Normal Transformed Correlations (Row & Column)" << endl 
       << Z << endl;
  */
  mcorr::mcorr_getPVals();
  k = mcorr::mcorr_BHY();
  /* cout << "-----------------------------------------------" << endl;
  cout << "Corresponding Sorted P-values" << endl << sortedPVals << endl;
  */
  cout << "-----------------------------------------------" << endl;
  cout << "Test results using Benjamini/Hochberg/Yekutieli" << endl;
  if (corrType == 0)
    cout << "for Pearson Correlations:" << endl;
  else if (corrType == 1)
    cout << "for Spearman Correlations:" << endl;
  else
    cout << "for Kendall Correlations:" << endl;
  cout << "Alpha = " << alpha << endl << endl;

  if(k == 0){
    cout << "Reject none of the null hypotheses" << endl << endl << endl;
    return;
  }
  
  cout << "Reject the null hypotheses of nonzero correlation corresponding to" << endl;
  for(int i = 0; i < k; i++){
    cout << "P-value_(" << i + 1 << "):  Vector " << (int)(sortedPVals[i][3] + 1) 
	 << " and Vector " << (int)(sortedPVals[i][2] + 1) << endl;
  }
  cout << endl << endl;
  
}


