#ifndef EXPERIMENT_HPP
#define EXPERIMENT_HPP

extern const int FLAG_TEST_CORR;
extern const int FLAG_TEST_MULT;
extern const int FLAG_TEST_PORT;
extern const int FLAG_TEST_ALL;

void runCorrelationTests(const char* filesWithStatus[], const int numberOfStatus,
    int flagTestEnable = FLAG_TEST_ALL, double Alpha = 1.0e-4);

#endif // EXPERIMENT_HPP