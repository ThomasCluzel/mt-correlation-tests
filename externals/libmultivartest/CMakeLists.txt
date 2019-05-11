cmake_minimum_required(VERSION 3.0)

project(libmultivartest
        LANGUAGES C CXX
        DESCRIPTION "Ismay's multivariate statistical tests")

# Compilation of the library => libmultivartest.a
add_library(multivartest STATIC
                        kendall2.c
                        mcorr.cpp mcorr.h
                        mmult.cpp mmult.h
                        mport.cpp mport.h
                        prob/prob.cpp prob/prob.hpp)
target_include_directories(multivartest SYSTEM PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
