#ifndef ALLEXPERIMENTS_HPP
#define ALLEXPERIMENTS_HPP

/**
 * List of (n, p) possible:
 * (16384, 2), (8192, 4), (4096, 8), (2048, 16), (1024, 32), (512, 64), (256, 128)
 * 2 sequences of 16384 numbers, 4 sequences of 8192 numbers...
 */

#include "experiment.hpp"

/**
 * \fn void runExperiments(const char* filesWithStatus[], int flagTestsEnable = FLAG_TEST_ALL, double alpha = 1.0e-4);
 * \brief launch the tests for all (n, p) possible
 * \param filesWithStatus array of names of files containing a status of Mersenne Twister
 * \param flagTestsEnable flag to select tests to run
 * \param alpha authorised margin of error of the tests
 */
void runExperiments(const char* filesWithStatus[], int flagTestsEnable = FLAG_TEST_ALL, double alpha = 1.0e-4);

/**
 * \fn void runControlExperiments();
 * \brief run all the control experiments on each (n, p) possible
 */
void runControlExperiments();

/**
 * \fn void runExperimentsStatusSeparatedByOneTrillionNumbers();
 * \brief run all the experiments on each (n, p) possible for the 128 first
 *        status in the mersenne_twister_states folder
 */
void runExperimentsStatusSeparatedByOneTrillionNumbers();

#endif // ALLEXPERIMENTS_HPP