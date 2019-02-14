// C. Ismay's includes
#include "libmultivartest/tntjama/tnt.h"
#include "libmultivartest/mcorr.h"
#include "libmultivartest/mmult.h"
#include "libmultivartest/mport.h"

// Include Mersenne Twister
extern "C" {
#include "libmt/mt19937ar.h"
}

#include "experiment.hpp"

// const flags
const int FLAG_TEST_CORR = 0b001;
const int FLAG_TEST_MULT = 0b010;
const int FLAG_TEST_PORT = 0b100;
const int FLAG_TEST_ALL  = 0b111;

// Constants
const int NP = 32768;

// Functions
void runCorrelationTests(const char* filesWithStatus[], const int numberOfStatus,
                         int flagTestEnable, int lagOrder, double alpha)
{
    const int n = NP/numberOfStatus, p = numberOfStatus;
    TNT::Array2D <long double> mat(n, p);

    /*** Initialize the sequences ***/
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

    /*** run the tests ***/

    //To avoid scientific notation
    //cout.setf(std::ios::fixed);

    if(flagTestEnable & FLAG_TEST_CORR)
    {
        mcorr mcorr1 = mcorr(n, p, mat, alpha);

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
        mmult mmult1 = mmult(n, p, mat, alpha);
        mmult1.mmult_LRT();
    }
    
    if(flagTestEnable & FLAG_TEST_PORT)
    {
        //Portmanteau tests for white noise
        if(lagOrder <= 0 || lagOrder > p)
            lagOrder = p; // default value of lagOrder is p
        mport mport1 = mport(n, p, mat, alpha);
        mport1.mport_centerMat();
        mport1.mport_portmanteauTests(lagOrder);
        mport1.mport_mahdiMcLeod(lagOrder);
    }
}