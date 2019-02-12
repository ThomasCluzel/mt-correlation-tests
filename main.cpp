#include <iostream>
#include "experiment.hpp"

// main program that runs tests
int main(int, char**)
{
    const char* statusFiles[] = { "mersenne_twister_states/mts000M000000", // 8 times the same status
                                  "mersenne_twister_states/mts000M000000", // should show a correlation
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
