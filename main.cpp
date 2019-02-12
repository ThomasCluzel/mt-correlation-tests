#include <iostream>
#include "allexperiments.hpp"

int main(int, char**)
{
    runControlExperiments();
    runExperimentsStatusSeparatedByOneTrillionNumbers();
    // TODO: run other experiments
    // - the second half (128 to 255)
    // - separate by 2 trillions (0 to 256)
    // - change the lag order to p/2
    return 0;
}
