#include <iostream>
#include "allexperiments.hpp"

using std::cout;

int main(int, char**)
{
    // control experiment
    cout << "Beginning of the control experiments (all status are the same)\n";
    runControlExperiments();

    // experiment
    cout << "\nBeginning of the experiments with status separated by one trillion numbers\n";
    runAllExperiments();

    // test the influence of distance between status
    cout << "\nBeginning of the experiments with status separated by two trillion numbers\n";
    runAllExperiments(2);

    // test the influence of the lag order
    cout << "\nBeginning of the experiments with a lag order of p/2\n";
    runExperimentsOnLagOrder();
    
    return 0;
}
