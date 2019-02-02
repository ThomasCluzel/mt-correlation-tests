cmake_minimum_required(VERSION 3.0)

project(libmultivartest LANGUAGES C CXX)

# Variables containing the sources files
file(GLOB SRC *.c *.cpp)
file(GLOB HEADER *.h)
file(GLOB PROB prob/*)
file(GLOB TNTJAMA tntjama/*)

# Compilation of the library (static) => libmultivartest
add_library(multivartest ${SRC} ${HEADER} ${PROB} ${TNTJAMA})
