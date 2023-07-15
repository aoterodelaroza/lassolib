#ifndef _LIBLASSO_H
#define _LIBLASSO_H

// Interface to the LIBLASSO. These interfaces correspond to the functions
// defined in c_functions.f90.

extern "C"{
void lasso_c(int rows, int cols, double *x, double *y, double t, double *w);
}

#endif
