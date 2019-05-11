cmake_minimum_required(VERSION 3.0)

project(libmt19937ar
        LANGUAGES C
        DESCRIPTION "Mersenne Twister PRNG by Makoto Matsumoto and Takuji Nishimura")

# Compilation of the library => libmt19937ar.a
add_library(mt19937ar STATIC mt19937ar.c mt19937ar.h)
target_include_directories(mt19937ar SYSTEM PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
