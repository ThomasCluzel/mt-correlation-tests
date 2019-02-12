#ifndef EXPERIMENT_HPP
#define EXPERIMENT_HPP

/**
 * list of flags usable by the `flagTestEnable` argument
 * of the function runCorrelationTests
 */
extern const int FLAG_TEST_CORR; // !< Enable the mcorr tests
extern const int FLAG_TEST_MULT; // !< Enable the mmult tests
extern const int FLAG_TEST_PORT; // !< Enable the mport tests
extern const int FLAG_TEST_ALL;  // !< Enable the all the tests

/**
 * \fn void runCorrelationTests(const char* filesWithStatus[], const int numberOfStatus,
            int flagTestEnable = FLAG_TEST_ALL, double Alpha = 1.0e-4);
 * \brief perform the statistical tests on a bunch of initial status
 * \param filesWithStatus the set of initial status (the names of the files)
 *        containing them
 * \param numberOfStatus size of the array of status
 * \param flagTestEnable flags to choose which test to enable
 * \param Alpha the maximum error authorised (results are given \f$ \pm Alpha \f$)
 */
void runCorrelationTests(const char* filesWithStatus[], const int numberOfStatus,
    int flagTestEnable = FLAG_TEST_ALL, double Alpha = 1.0e-4);

#endif // EXPERIMENT_HPP