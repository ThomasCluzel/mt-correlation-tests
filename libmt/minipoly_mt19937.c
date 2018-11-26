#include <NTL/GF2X.h>
#include <NTL/vec_GF2.h>
#include <cstdlib>
#include <fstream>
#include "mt19937ar.h"
#include "shared.h"

#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

#define MEXP 19937 /* the dimension of the state space */

/* computes the minimal polynomial of the linear recurrence */
void comp_mini_poly (void)
{
  int i;
  vec_GF2 v (INIT_SIZE, 2 * MEXP);

  for (i=0; i<2*MEXP; i++)
    v[i] = genrand_int32() & 0x01ul ;

  MinPolySeq (phi, v, MEXP);
}

/* computes the t^J mod phi(t) */
void comp_jump_rem (long jump_step)
{
  /* changed by saito 2013.1.25 */
  //GF2X f;
  //  SetCoeff (f, jump_step, 1);
  //  g = f % phi;
  PowerXMod(g, jump_step, phi);
  /* changed by saito 2013.1.25 */
}

/*
int main (void)
{
  int i, a=0;
  long jump_step = 1234567;  the number of steps of jumping ahead
  unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
  ofstream fout;

  init_by_array(init, length);

  comp_mini_poly ();
  comp_jump_rem (jump_step);

  fout.open ("clist_mt19937.txt", ios::out);
  if (!fout)
    return -1;

  for (i=MEXP-1; i>-1; i--)
    fout << coeff (g, i);

  fout.close();

  return 0;
}*/

