cmake_minimum_required(VERSION 3.0)

# These are the tests
project(tests)

# Include the library
link_directories(../libmt)

# The test programs and library link
add_executable(test_mt mtTest.c)
target_link_libraries(test_mt m)
target_link_libraries(test_mt mt19937ar)
add_executable(test_mt_status mtStateTest.c)
target_link_libraries(test_mt_status mt19937ar)

# CTest
enable_testing()

# The tests
add_test(test1 test_mt)
add_test(test2 test_mt_status) # this test is very long
