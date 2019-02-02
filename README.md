# Mersenne Twister: Correlation of Sub-sequences

## Motivation

Is the Mersenne Twister Random Number Generator a good RNG for
multi-threaded applications?

The goal of this project is to determinate if one may use Mersenne
Twister in order to generate non-correlated sequences of random numbers.

## How do we determinate this?

We test sub-sequences generated by Mersenne Twister with inter-streams correlation
tests of C. Ismay (described in his PhD thesis "_Testing Independence of Parallel
Pseudorandom Number Streams Incorporating the Data’s Multivariate Nature_").

## Description of the repository

This repository contains:
- a libmt folder: it contains the Mersenne Twister PRNG as a library
- a libmultivartest folder: it contains C. Ismay's correlation tests
- a tests folder: it contains unit tests that validate the implementation of Mersenne Twister
- a mersenne_twister_states folder: it contains status of the PRNG separated by 1000 billion numbers.
- files that use the tests to test the existence of correlations in sequences of Mersenne Twister

### How to compile

This project uses CMake to offer the largest possible compatibility.
1. To install CMake have a look at this page
[https://cmake.org/download/](https://cmake.org/download/)
1. Set current directory to the root of this project(`cd path/to/folder/`)
1. Then run `cmake -G "Unix Makefiles"` to use linux make or
`cmake -G "MinGW Makefiles"` to use MinGW make (if you want to use an
IDE look for its generator in the list of generators)
1. Compile the project either with the tool you are used to or by running
`cmake --build .`
1. Finally you can run the binary

### How to run the unit tests

After having compiled the project, go in the "tests" folder
(`cd tests`) and run the command `ctest -V`.

Take note that you may want to disable the unit test that
check the gap between two status because it is very long.
However you should not trust us and run the test to verify that
the status are separated by enough numbers.

### How to run the statistical tests

Run the program `./prog` (linux) or `prog.exe` (windows).

### "What if I don't want to install cmake?"

"No problem", for an easier reproducibility, you can use docker.
However, if you don't know anything about docker or virtualisation, you
may want to compile the project yourself rather than using docker.

To compile and run the code with docker:
- Compile the code (and run unit tests) by building the image:
`docker build -t mt-correlation-tests .`
- Run the statistical tests by running a container:
`docker run --rm mt-correlation-tests`

:warning: It is mandatory to have an internet connection the first time you
build the container because it needs to download the base image namely
the one containing gcc and cmake.

-------------------------------------------------------------------------

Bruno JOUSSE & Thomas CLUZEL
