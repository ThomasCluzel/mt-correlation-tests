#include "experiment.hpp"

void runControlExperiments()
{
    // initialize the status
    const char* initialStatus = "mersenne_twister_states/mts000M000000";
    const char* statusFiles[128];
    for(int i=0; i<128; i++)
        statusFiles[i] = initialStatus;

    // run the control experiments
    for(int numberOfStatus=2; numberOfStatus!=128; numberOfStatus*=2)    
    {
        runCorrelationTests(statusFiles, numberOfStatus, FLAG_TEST_CORR | FLAG_TEST_MULT);
    }
}