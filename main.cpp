// std C++ includes
#include <iostream>
#include <fstream>

// C. Ismay's includes
#include "libmultivartest/tntjama/tnt.h"
#include "libmultivartest/mcorr.h"
#include "libmultivartest/mmult.h"
#include "libmultivartest/mport.h"
#include "libmultivartest/prob/prob.cpp"

// Include Mersenne Twister
extern "C" {
#include "libmt/mt19937ar.h"
}

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::flush;
using std::string;
using std::exit;

// prototypes of functions
void runCorrelationTests(const char* filesWithStatus[], const int numberOfStatus);

// main program that runs tests
int main(int, char**)
{
    const char* statusFiles[] = { "mersenne_twister_states/mts000M000000", // 10 times the same status
                                  "mersenne_twister_states/mts000M000000", // should show a correlation
                                  "mersenne_twister_states/mts000M000000",
                                  "mersenne_twister_states/mts000M000000",
                                  "mersenne_twister_states/mts000M000000",
                                  "mersenne_twister_states/mts000M000000",
                                  "mersenne_twister_states/mts000M000000",
                                  "mersenne_twister_states/mts000M000000",
                                  "mersenne_twister_states/mts000M000000",
                                  "mersenne_twister_states/mts000M000000"};
    const int numberOfStatus = sizeof(statusFiles) / sizeof(statusFiles[0]);

    // control experiment
    runCorrelationTests(statusFiles, numberOfStatus);

    return 0;
}

/**
 * \fn void runCorrelationTests(string filesWithStatus[])
 * \brief A function to perform the statistical
 *        tests on a bunch of initial status
 * \param filesWithStatus an array of file names that contain status of MT
 */
void runCorrelationTests(const char* filesWithStatus[], const int numberOfStatus)
{
    // first open an empty file to store values generated in it
    const char* tmpfilename = "tmp.txt";
    ofstream fileContainingNumbers(tmpfilename, ios_base::out | ios_base::trunc);
    if(!fileContainingNumbers.is_open())
    {
        cerr << "Error: could not open file " << tmpfilename << endl;
        exit(1);
    }

    // set all the properties required to perform the test
    const int n = 100000/numberOfStatus, p = numberOfStatus; // original values were n=10000, p=10
    TNT::Array2D <long double> mat(n, p);
    double Alpha;

    for(int i=0; i<p; i++) // for each status
    {
        // restore the status into MT status
        restoreStatus(filesWithStatus[i]);
        for(int j=0; j<n; j++)
        {
            // fill the file with a sequence generated from this status
            fileContainingNumbers << genrand_real1() << "\n";
        }
    }
    fileContainingNumbers << flush;
    fileContainingNumbers.close();

    // Generate the matrix from the file
    ifstream fin;
    fin.open(tmpfilename);
    if(!fin.is_open()){
        cout<< "Error opening input file!" << endl;
        exit(2);
    }
    //Input file one row after another
    for(int i = 0; i < n; i++)
        for(int j = 0; j < p; j++)
            fin >> mat[i][j];
    fin.close();
    remove(tmpfilename);

    // run all the tests

    //std::cout << "Enter alpha value: " << endl;
    //std::cin >> Alpha;
    Alpha = 0.01;

    mcorr mcorr1 = mcorr(n, p, mat, Alpha);

    //To avoid scientific notation
    //cout.setf(std::ios::fixed);

    //Pearson
    mcorr1.mcorr_pairCorr(0);
    //Spearman
    mcorr1.mcorr_pairCorr(1);
    //Kendall
    mcorr1.mcorr_pairCorr(2);

    //Testing if correlation matrix = identity
    mmult mmult1 = mmult(n, p, mat, Alpha);
    mmult1.mmult_LRT();

    //Portmanteau tests for white noise with lag = p
    mport mport1 = mport(n, p, mat, Alpha);
    mport1.mport_centerMat();
    mport1.mport_portmanteauTests(p);
    mport1.mport_mahdiMcLeod(p);
}