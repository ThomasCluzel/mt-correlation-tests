#include <stdio.h>
#include <stdlib.h>
#include "libmt/shared.h"
#include "libmt/jump_mt19937.h"
#include "libmt/mt19937ar.h"

int main(int argc, char* argv[])
{
    unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
    init_by_array(init, length);
    
    int i, j, jumpLength, deg, c;
    unsigned long *pf;
    State *ss1, *ss2, ss3;
    FILE *fin;

    if(argc<2){
        puts("Vous devez saisir le nombre d'itérations en paramètre.");
        return 1;
    }
    jumpLength=strtol(argv[1], NULL, 10);

	//Calculate polynomial
    comp_mini_poly ();
	comp_jump_rem (jumpLength);
	
    pf = (unsigned long *)calloc(P_SIZE, sizeof(unsigned long));
    
    //Get polynomial coefficient and give it to pf
    for (i=MEXP-1; i>-1; i--){
		c = coeff(g, i);
		set_coef(pf, i, c);
    }
    
    //Set initial state
    for(i=0; i<624; i++){
		s0->s[i] = mt[i];
	}
	
	//Jump manually
    for(i=0; i<jumpLenght; i++){
        printf("%lf \n", genrand_real1()) ;
    }
    
    /*Computes jumping ahead with standard Horner method */ 
    ss1 = horner1(pf, &s0);
    
    /*Computes jumping ahead with Sliding window algorithm*/
    gen_vec_h(&s0);
    ss2 = calc_state(pf, &s0);
    
    puts("Comparaison de l'état:");
    puts("mot |   MT   | Horner | Sliding");

    for(i=0; i<624; i++){
        printf("%d | %ld | %ld | %ld\n", i, mt[i]; ss1.s[i], ss2.s[i]);
    }

    return 0;
