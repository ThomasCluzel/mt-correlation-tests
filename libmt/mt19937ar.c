/* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          
   Copyright (C) 2005, Mutsuo Saito,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

#include <NTL/GF2X.h>
#include <NTL/vec_GF2.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

NTL_CLIENT

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

/* Parameter for computing the minimal polynomial */
#define MEXP 19937 /* the dimension of the state space */

/* parameters for computing Jump */
#define W_SIZE 32 /* size of unsigned long */
//#define MEXP 19937
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

/* initial state */
//static State s0;

static unsigned long h[LL]; 
static State vec_h[LL];

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */


/**
 * File minipoly_mt
 */

GF2X phi; /* phi is the minimal polynomial */
GF2X g; /* g(t) is used to store t^J mod phi(t) */

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

int compute_clist()
{
  int i, a=0;
  long jump_step = 1234567; /* the number of steps of jumping ahead */
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
}

/**
 * File jump_mt
 */
int main_old(void)
{
  unsigned long *pf;
  State *ss1, *ss2, ss3;
  int i, j, deg;
  char c;
  FILE *fin;
  State s0;

  pf = (unsigned long *)calloc(P_SIZE, sizeof(unsigned long));

  /* read the file clist.txt, and set the coefficient */
  if ((fin = fopen("clist_mt19937.txt","r")) == NULL){
    printf("File read error.\n");
    exit(1);
  }

  for (i=MEXP-1; i>-1; i--){
    c = fgetc(fin);
    if (c == '1')
      set_coef(pf, i, 1);
  }

  /* computes jumping ahead with standard Horner method */    
  ss1 = horner1(pf, &s0);
    
  /* computes jumping ahead with Sliding window algorithm */
  gen_vec_h(&s0);
  ss2 = calc_state(pf, &s0);

  if (compare_state(ss1, ss2) != 0)
    printf("error\n");
  
  return(0);  
}

/* generate random polynomial */
void random_poly(unsigned long *p)
{
  int i;

  for (i=0; i<MEXP; i++)
      set_coef(p, i, rand() % 2);
}

/* copy state ss */
void copy_state(State *ts, State *ss)
{
  int i;

  for(i=0; i<N; i++)
    ts->s[i] = ss->s[i];

  ts->ptr = ss->ptr;
}

/* compute pf(ss) using standard Horner method */
State *horner1(unsigned long *pf, State *ss)
{
  int i=MEXP-1;
  State *temp;

  temp = (State *)calloc(1, sizeof(State));

  while(get_coef(pf,i) == 0)
    i--;

  if (i > 0){
    copy_state(temp, ss);
    gen_next(temp);
    i--;
    for (; i>0; i--){
      if (get_coef(pf,i) != 0)
	add_state(temp, ss);
      else
	;
      gen_next(temp);
    }
    if (get_coef(pf,0) != 0)
      add_state(temp, ss);
    else
      ;
  }
  else if(i == 0)
    copy_state(temp,ss);
  else 
    ;
  
  return(temp);

}

/* compute next state */

/* add s2 to s1 */
void add_state(State *s1, State *s2) 
{
  int i, pt1=s1->ptr, pt2=s2->ptr;
  
  if (pt2-pt1 >= 0){
    for (i=0; i<N-pt2; i++)
      s1->s[i+pt1] ^= s2->s[i+pt2];
    for(; i<N-pt1; i++)
      s1->s[i+pt1] ^= s2->s[i+(pt2-N)];
    for(; i<N; i++)
      s1->s[i+(pt1-N)] ^= s2->s[i+(pt2-N)];
  }
  else{
    for (i=0; i<N-pt1; i++)
      s1->s[i+pt1] ^= s2->s[i+pt2];
    for(; i<N-pt2; i++)
      s1->s[i+(pt1-N)] ^= s2->s[i+pt2];
    for(; i<N; i++)
      s1->s[i+(pt1-N)] ^= s2->s[i+(pt2-N)];
  }
}

/* compute pf(ss) using Sliding window algorithm */
State *calc_state(unsigned long *pf, State *ss)
{
  State *temp1;
  int i=MEXP-1, j, digit, skip=0;

  temp1 = (State *)calloc(1, sizeof(State));

  while (get_coef(pf,i) == 0)
    i--;

  for (; i>=QQ; i--){
    if (get_coef(pf,i) != 0){
      for (j=0; j<QQ+1; j++)
	gen_next(temp1);
      digit = 0; 
      for (j=0; j<QQ; j++)
	digit = (digit << 1) ^ get_coef(pf, i-j-1);
      add_state(temp1, &vec_h[digit]);
      i -= QQ;
    }
    else
      gen_next(temp1);
  }
  
  for (; i>-1; i--){
    gen_next(temp1);
    if (get_coef(pf,i) == 1)
      add_state(temp1, ss);
    else
      ;
  }
  return(temp1);
}

/* compare s1 with s2 */
/* if s1 is different from s2, then return 1 , else return 0 */
int compare_state(State *s1, State *s2)
{
  int i, d = 0;

  if ((s1->s[s1->ptr] & 0x80000000ul) != (s2->s[s2->ptr] & 0x80000000ul))
    d = 1;

  for (i=1; i<N; i++)
    if (s1->s[(i+s1->ptr)%N] != s2->s[(i+s2->ptr)%N]){
      d = 1;
      break;
    }

  return(d);
}

/* generate Gray code */
void gray_code(void)
{
  unsigned int i, j=1, l=1, term = LL;

  h[0] = 0;
  
  for (i=1; i<=QQ; i++){
      l = (l << 1);
      term = (term >> 1);
      for (; j<l; j++)
	  h[j] = h[l-j-1] ^ term; 
  }
}

/* compute h(f)ss where h(t) are exact q-degree polynomials, */
/* f is the tarnsition function, and ss is the initial state */
/* the results are stored in vec_h[0] , ... , vec_h[LL-1]    */
void gen_vec_h(State *ss)
{
  int i;
  unsigned long k, g;
  State v;

  gray_code();

  copy_state(&vec_h[0], ss);

  for(i=0; i<QQ; i++)
    gen_next(&vec_h[0]);

  for (i=1; i<LL; i++){
    copy_state(&v, ss);
    g = h[i] ^ h[i-1];
    for (k=1; k<g; k=(k<<1))
      gen_next(&v);
    copy_state(&vec_h[h[i]], &vec_h[h[i-1]]);
    add_state(&vec_h[h[i]], &v);
  }
}


/* 32-bits function */ 
/* return the i-th coefficient of the polynomial pf */
unsigned long get_coef(unsigned long *pf, unsigned int deg)
{
  if ((pf[deg >> 5] & (LSB << (deg & 0x1ful))) != 0) 
    return(1);
  else
    return(0);
}

/* 32-bit function */
/* set the coefficient of the polynomial pf with v */
void set_coef(unsigned long *pf, unsigned int deg, unsigned long v)
{
  if (v != 0)
    pf[deg >> 5] ^= (LSB << (deg & 0x1ful));
  else
    ;
}

/* next state generating function */
void gen_next(State *ss)
{
  int num;
  unsigned long y;
  static unsigned long mag02[2]={0x0ul, MATRIX_A};

  num = ss->ptr;
  if (num < N-M){
      y = (ss->s[num] & UPPER_MASK) | (ss->s[num+1] & LOWER_MASK);
      ss->s[num] = ss->s[num+M] ^ (y >> 1) ^ mag02[y % 2];
      ss->ptr++;
  }
  else if (num < N-1){
      y = (ss->s[num] & UPPER_MASK) | (ss->s[num+1] & LOWER_MASK);
      ss->s[num] = ss->s[num+(M-N)] ^ (y >> 1) ^ mag02[y % 2];
      ss->ptr++;
  }
  else if (num == N-1){
      y = (ss->s[N-1] & UPPER_MASK) | (ss->s[0] & LOWER_MASK);
      ss->s[N-1] = ss->s[M-1] ^ (y >> 1) ^ mag02[y % 2];
      ss->ptr = 0;
  }
}


/**
 * Our functions
 */
void getState(State* s)
{
  memcpy(s->s, mt, N*sizeof(unsigned long));
  s->ptr = mti;
}
void setState(State* s)
{
  memcpy(mt, s->s, N*sizeof(unsigned long));
  mti = s->ptr;
}
void jump_ahead(long steps)
{
  /*
  State s0;
  State* ss2;
  std::stringstream buf;
  std::string buf2;
  //char c;
  int i, j;
  unsigned long pf[P_SIZE];
  //pf = (unsigned long *)calloc(P_SIZE, sizeof(unsigned long));
  for(i=0 ; i<P_SIZE ; i++)
    pf[i] = 0;

  comp_mini_poly();
  comp_jump_rem(steps);

  for (i=MEXP-1; i>-1; i--)
    buf << coeff (g, i);
  buf2 = buf.str();
  for (i=MEXP-1, j=0; i>-1; i--, j++)
  {
    //c = fgetc(fin);
    //c = buf2[j];
    std::cout << buf2 << std::endl;
    if (buf2[j] == '1')
      set_coef(pf, i, 1);
  }

  getState(&s0);
  gen_vec_h(&s0);
  ss2 = calc_state(pf, &s0);
  setState(ss2);
  free(ss2);
  */
  int i, a=0;
  long jump_step = steps; /* the number of steps of jumping ahead */
  //unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
  ofstream fout;

  //init_by_array(init, length);

  comp_mini_poly ();
  comp_jump_rem (jump_step);

  fout.open ("clist_mt19937.txt", ios::out);
  if (!fout)
    exit(-1);

  for (i=MEXP-1; i>-1; i--)
    fout << coeff (g, i);

  fout.close();

  // then read the file

  unsigned long *pf;
  State *ss1, *ss2, ss3;
  State s0;
  int j, deg;
  char c;
  FILE *fin;

  pf = (unsigned long *)calloc(P_SIZE, sizeof(unsigned long));

  /* read the file clist.txt, and set the coefficient */
  if ((fin = fopen("clist_mt19937.txt","r")) == NULL){
    printf("File read error.\n");
    exit(1);
  }

  for (i=MEXP-1; i>-1; i--){
    c = fgetc(fin);
    if (c == '1')
      set_coef(pf, i, 1);
  }
  fclose(fin);

  getState(&s0);

  /* computes jumping ahead with standard Horner method */    
  ss1 = horner1(pf, &s0);
    
  /* computes jumping ahead with Sliding window algorithm */
  gen_vec_h(&s0);
  ss2 = calc_state(pf, &s0);

  if (compare_state(ss1, ss2) != 0)
    printf("error the states are different\n");

  setState(ss2);

  free(ss1);
  free(ss2);
  free(pf);

  remove("clist_mt19937.txt");
}

// just a quick test
int main()
{
  long i, jump_step = 10; /* the number of steps of jumping ahead */
  unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
  State si;
  unsigned long output[10], o, nbGen=10;

  init_by_array(init, length);

  getState(&si); // save the initial state
  printf("mti=%d\n", mti);

  for(i=0 ; i<jump_step ; i++) // skip x numbers
    genrand_int32();
  printf("mti=%d\n", mti);

  for(i=0 ; i<nbGen ; i++)
    output[i] = genrand_int32();

  setState(&si); // reset the state
  printf("mti=%d\n", mti);

  jump_ahead(jump_step); //jump x steps
  printf("mti=%d\n", mti);

  for(i=0 ; i<nbGen ; i++)
  {
    o = genrand_int32();
    if(o != output[i])
      printf("Error: expected %lu, got %lu\n", output[i], o);
  }
  printf("End \n");

  return 0;
}

/**
 * TODO:
 * - find why it failed
 * - move the main function in a file in the tests folder
 * - update readme
 * - update CMakeLists.txt
 */