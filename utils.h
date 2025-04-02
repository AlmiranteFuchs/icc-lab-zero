#ifndef __UTILS_H__
#define __UTILS_H__

#include "DoubleType.h"
#include <stdlib.h>
#include <time.h>


typedef double real_t;
typedef double rtime_t;
typedef char * string_t;

double timestamp(void);

real_t abs_value(real_t x);

int64_t ulp_distance(Double_t x, Double_t y);

#endif // __UTILS_H__

