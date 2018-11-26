#ifndef JUMP_MT19937_H
#define JUMP_MT19937_H

#include <stdio.h>
#include <stdlib.h>
#include "shared.h"

/* parameters of Mersenne Twister */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

/* parameters for computing Jump */
#define W_SIZE 32 /* size of unsigned long */
#define MEXP 19937
#define P_SIZE ((MEXP/W_SIZE)+1)
#define LSB 0x00000001ul
#define QQ 7
#define LL 128  /* LL = 2^(QQ) */

typedef struct{
  unsigned long s[N];
  int ptr;
}State;

/* operations for polynomials */
void random_poly(unsigned long *);
void gray_code(void);
unsigned long get_coef(unsigned long *, unsigned int);
void set_coef(unsigned long *, unsigned int, unsigned long);

/* operations for states */
State *calc_state(unsigned long *, State *);
void copy_state(State *, State *);
State *horner1(unsigned long *, State *);
void gen_next(State *);
void add_state(State *, State *);
int compare_state(State *, State *);
void gen_vec_h(State *);

#endif
