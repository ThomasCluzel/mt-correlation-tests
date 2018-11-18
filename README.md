# Mersenne Twister Correlation of Subsequences

## Motivation

Is the Mersenne Twister Random Number Generator a good RNG for
multi-threaded application?
The goal of this project is to determinate how one can use the
MT RNG in order to generate non-correlated sequences of random number.

## How

We initialize MT.
We use the jump ahead algorithm to jump from one state to another.
By doing this we generate sequences of random numbers extracted from
the (huge) period of MT.
We test these sequences with the TestU01 inter-streams correlation tests.

-------------------------------------------------------------------------

Bruno JOUSSE & Thomas CLUZEL
