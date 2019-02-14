#include <iostream>
#include "allexperiments.hpp"

int main(int, char**)
{
    runControlExperiments();
    runAllExperiments();
    runAllExperiments(2); // influence of distance between status
    // TODO: run other experiments
    // - change the lag order to p/2
    return 0;
}
