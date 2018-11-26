#ifndef MINIPOLY_MT19937__H
#define MINIPOLY_MT19937__H

#include <NTL/GF2X.h>
#include <NTL/vec_GF2.h>
#include <cstdlib>
#include <fstream>
#include "mt19937ar.h"
#include "shared.h"

NTL_CLIENT

void comp_mini_poly (void);
void comp_jump_rem (long jump_step);

#endif
