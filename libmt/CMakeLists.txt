cmake_minimum_required(VERSION 3.0)

project(libmt19937ar)

set(SRC mt19937ar.c)
set(HEADER mt19937ar.h)

# Compilation of the library (static) => libmt19937ar
add_library(mt19937ar ${SRC} ${HEADER})
