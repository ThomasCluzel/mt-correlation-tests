// std C++ includes
#include <iostream>

// C. Ismay's includes
#include "libmultivartest/tntjama/tnt.h"
#include "libmultivartest/mcorr.h"
#include "libmultivartest/mmult.h"
#include "libmultivartest/mport.h"
#include "libmultivartest/prob/prob.cpp" // weird but it is in Ismay's thesis

// Include Mersenne Twister
extern "C" {
#include "libmt/mt19937ar.h"
}

using std::cout;
using std::endl;
using std::exit;

// constants
const int FLAG_TEST_CORR = 0b001;
const int FLAG_TEST_MULT = 0b010;
const int FLAG_TEST_PORT = 0b100;
const int FLAG_TEST_ALL  = 0b111;

// prototypes of functions
void runCorrelationTests(const char* filesWithStatus[], const int numberOfStatus, int flagTestEnable = FLAG_TEST_ALL);


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
                                  "mersenne_twister_states/mts000M000000" };
    const int numberOfStatus = sizeof(statusFiles) / sizeof(statusFiles[0]);

    // control experiment
    runCorrelationTests(statusFiles, numberOfStatus, FLAG_TEST_CORR | FLAG_TEST_MULT);

    return 0;
}

/* A function to perform the statistical */
/* tests on a bunch of initial status    */
void runCorrelationTests(const char* filesWithStatus[], const int numberOfStatus, int flagTestEnable)
{
    const int n = 100000/numberOfStatus, p = numberOfStatus; // original values were n=10000, p=10
    TNT::Array2D <long double> mat(n, p);
    double Alpha;

    for(int s=0; s<p; s++) // for each status
    {
        // set MT status to this status
        restoreStatus(filesWithStatus[s]);
        for(int r=0; r<n; r++) // for each random number in the sequence
        {
            // populate a column of the matrix with the random numbers 
            mat[r][s] = genrand_real1();
        }
    }

    /*** run all the tests ***/

    //To avoid scientific notation
    //cout.setf(std::ios::fixed);

    //std::cout << "Enter alpha value: " << endl;
    //std::cin >> Alpha;
    Alpha = 0.01;

    if(flagTestEnable & FLAG_TEST_CORR)
    {
        mcorr mcorr1 = mcorr(n, p, mat, Alpha);

        //Pearson
        mcorr1.mcorr_pairCorr(0);
        //Spearman
        mcorr1.mcorr_pairCorr(1);
        //Kendall
        mcorr1.mcorr_pairCorr(2);
    }
    
    if(flagTestEnable & FLAG_TEST_MULT)
    {
        //Testing if correlation matrix = identity
        mmult mmult1 = mmult(n, p, mat, Alpha);
        mmult1.mmult_LRT();
    }
    
    if(flagTestEnable & FLAG_TEST_PORT)
    {
        //Portmanteau tests for white noise with lag = p
        mport mport1 = mport(n, p, mat, Alpha);
        mport1.mport_centerMat();
        mport1.mport_portmanteauTests(p);
        mport1.mport_mahdiMcLeod(p);
    }
}