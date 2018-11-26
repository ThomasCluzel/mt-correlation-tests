#include "jump_ahead.h"

void jump_ahead_comparison(int jumpLength){
	int i, j, c;
    unsigned long *pf;
    State *ss1, *ss2, ss3;

	//Calculate polynomial
    comp_mini_poly ();
	comp_jump_rem (jumpLength);
	
    pf = (unsigned long *)calloc(P_SIZE, sizeof(unsigned long));
    
    //Get polynomial coefficient and give it to pf
    for (i=MEXP-1; i>-1; i--){
		c = coeff(g, i);
		set_coef(pf, i, c);
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
	
}

unsigned long* jump_ahead_manually(int jumpLength){
	for(i=0; i<jumpLenght; i++){
        //printf("%lf \n", genrand_real1()) ;
    }
    
    return mt;
}

unsigned long* jump_ahead_horner(int jumpLength){
	int i, c;
    unsigned long *pf;
    State *ss;

	//Calculate polynomial
    comp_mini_poly();
	comp_jump_rem (jumpLength);
	
    pf = (unsigned long *)calloc(P_SIZE, sizeof(unsigned long));
    
    //Get polynomial coefficient and give it to pf
    for (i=MEXP-1; i>-1; i--){
		c = coeff(g, i);
		set_coef(pf, i, c);
    }
    
    /*Computes jumping ahead with standard Horner method */ 
    ss = horner1(pf, &s0);

	return ss->s;
}

unsigned long* jump_ahead_sliding(int jumpLength){
	int i, c;
    unsigned long *pf;
    State *ss;

	//Calculate polynomial
    comp_mini_poly();
	comp_jump_rem (jumpLength);
	
    pf = (unsigned long *)calloc(P_SIZE, sizeof(unsigned long));
    
    //Get polynomial coefficient and give it to pf
    for (i=MEXP-1; i>-1; i--){
		c = coeff(g, i);
		set_coef(pf, i, c);
    }
    
    /*Computes jumping ahead with standard Horner method */ 
    ss = calc_state(pf, &s0);

	return ss->s;
}
