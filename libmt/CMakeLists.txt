cmake_minimum_required(VERSION 3.0)

project(libmt19937ar LANGUAGES C)

# Variables containing the sources files
file(GLOB SRC *.c)
file(GLOB HEADER *.h)

# Compilation of the library (static) => libmt19937ar
add_library(mt19937ar ${SRC} ${HEADER})
