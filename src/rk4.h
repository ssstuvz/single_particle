#ifndef RK4_H
#define RK4_H
#include "common.h"

void derivatives_electron(my_float t, my_float y[], my_float dy_dt[]);

void rk4_step(my_float y_in[], int number_of_eqns, my_float t,
 my_float dt, my_float y_out[], void (*derivatives)(my_float, my_float[], my_float[]));

#endif