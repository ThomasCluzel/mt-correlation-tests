#include <iostream>
#include "allexperiments.hpp"

const int MAX_NUMBER_OF_STATUS = 128;

using namespace std;

void runExperiments(const char* filesWithStatus[MAX_NUMBER_OF_STATUS], int flagTestsEnable)
{
    for(int numberOfStatus=2; numberOfStatus<=MAX_NUMBER_OF_STATUS; numberOfStatus*=2) // for each (n, p) possible
    {
        cout << "\n\n\nExperiments with " << numberOfStatus << " status\n\n";
        runCorrelationTests(filesWithStatus, numberOfStatus, flagTestsEnable); // run tests
    }
}

void runControlExperiments()
{
    // initialize the status
    const char* initialStatus = "mersenne_twister_states/mts000M000000";
    const char* statusFiles[MAX_NUMBER_OF_STATUS];
    for(int i=0; i<MAX_NUMBER_OF_STATUS; i++)
        statusFiles[i] = initialStatus;

    // run the control experiments
    runExperiments(statusFiles, FLAG_TEST_CORR | FLAG_TEST_MULT);
}

void runAllExperiments(int step)
{
    // initialize the status
    const char* baseName = "mersenne_twister_states/mts000M000";
    char* statusFiles[MAX_NUMBER_OF_STATUS]; // array of pointers, not array of arrays
    for(int i=0; i<MAX_NUMBER_OF_STATUS; i++)
    {
        statusFiles[i] = (char*)malloc(38 * sizeof(char));
        if(statusFiles[i] == NULL)
        {
            cerr << "malloc error" << endl;
            exit(1);
        }
        sprintf(statusFiles[i], "%s%03d", baseName, i*step);
    }

    // run the all the experiments
    runExperiments((const char**)statusFiles);

    // free memory
    for(int i=0; i<MAX_NUMBER_OF_STATUS; i++)
        free(statusFiles[i]);
}

void runExperimentsOnLagOrder()
{
    // initialize the status
    const char* baseName = "mersenne_twister_states/mts000M000";
    char* statusFiles[MAX_NUMBER_OF_STATUS]; // array of pointers, not array of arrays
    for(int i=0; i<MAX_NUMBER_OF_STATUS; i++)
    {
        statusFiles[i] = (char*)malloc(38 * sizeof(char));
        if(statusFiles[i] == NULL)
        {
            cerr << "malloc error" << endl;
            exit(1);
        }
        sprintf(statusFiles[i], "%s%03d", baseName, i);
    }

    // run the experiments with lag order = numberOfStatus / 2
    for(int numberOfStatus=2; numberOfStatus<=MAX_NUMBER_OF_STATUS; numberOfStatus*=2) // for each (n, p) possible
    {
        cout << "\n\n\nExperiments with " << numberOfStatus << " status\n\n";
        runCorrelationTests((const char**)statusFiles, numberOfStatus, FLAG_TEST_ALL, numberOfStatus/2); // run tests
    }

    // free memory
    for(int i=0; i<MAX_NUMBER_OF_STATUS; i++)
        free(statusFiles[i]);
}