#ifndef JUMP_AHEAD__H
#define JUMP_AHEAD__H

#include "jump_mt19937.h"
#include "mt19937ar.h"
#include "minipoly_mt19937.h"
#include "shared.h"

NTL_CLIENT

void jump_ahead_comparison(int jumpLength);

unsigned long* jump_ahead_manually(int jumpLength);

unsigned long* jump_ahead_horner(int jumpLength);

unsigned long* jump_ahead_sliding(int jumpLength);

#endif
