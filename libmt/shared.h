#ifndef SHARED__H
#define SHARED_H

extern NTL_CLIENT

extern static unsigned long mt[N]; /* the array for the state vector  */
extern static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

extern GF2X phi; /* phi is the minimal polynomial */
extern GF2X g; /* g(t) is used to store t^J mod phi(t) */

/* initial state */
extern State s0;

extern unsigned long h[LL]; 
extern State vec_h[LL];

#endif
