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

#include <stdio.h>
#include <math.h>
#include "../libmt/mt19937ar.h"

#define FILENAME_CORRECT_OUTPUT "mt19937ar.out"
#define BUF_SIZE 512
#define EPSILON 0.00000001

int main(void)
{
   char buf[BUF_SIZE] = "";
   unsigned long random_integer = 0;
   unsigned long random_integer_correct = 0;
   double random_floating = 0.0;
   double random_floating_correct = 0.0;
   FILE *file_correct_output = NULL;

   int i;
   unsigned long init[4] = {0x123, 0x234, 0x345, 0x456}, length = 4;

   init_by_array(init, length);

   file_correct_output = fopen(FILENAME_CORRECT_OUTPUT, "r");
   if(file_correct_output == NULL)
   {
      fprintf(stderr, "Error: couldn't open file %s\n", FILENAME_CORRECT_OUTPUT);
      return 1; // test failed
   }

   //printf("1000 outputs of genrand_int32()\n");
   fgets(buf, BUF_SIZE, file_correct_output);
   for (i = 0; i < 1000; i++)
   {
      //printf("%10lu ", genrand_int32());
      random_integer = genrand_int32();
      fscanf(file_correct_output, "%lu", &random_integer_correct);
      if( random_integer != random_integer_correct )
      {
         fprintf(stderr, "Error: generated %lu expected %lu\n", random_integer, random_integer_correct);
         fclose(file_correct_output);
         return 1; // test failed
      }
      /*if (i % 5 == 4)
         printf("\n");*/
   }

   fgets(buf, BUF_SIZE, file_correct_output); // end of the line of integers
   fgets(buf, BUF_SIZE, file_correct_output); // the blank line between integers and reals
   fgets(buf, BUF_SIZE, file_correct_output); // the line introducing reals

   //printf("\n1000 outputs of genrand_real2()\n");
   for (i = 0; i < 1000; i++)
   {
      //printf("%10.8f ", genrand_real2());
      random_floating = genrand_real2();
      fscanf(file_correct_output, "%lf", &random_floating_correct);
      if( fabs(random_floating-random_floating_correct) > EPSILON )
      {
         fprintf(stderr, "Error: generated %10.8f expected %10.8f\n", random_floating, random_floating_correct);
         fclose(file_correct_output);
         return 1;
      }
      /*if (i % 5 == 4)
         printf("\n");*/
   }

   fclose(file_correct_output);

   return 0;
}
