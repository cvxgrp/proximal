#ifndef PROXIMAL_H_GUARD
#define PROXIMAL_H_GUARD

#include <math.h>
#include <float.h>

void project_box(double *x, const int n, const double *l, const double *u);
void project_pos(double *x, const int n);
void project_soc(double *x, const int n);
void prox_separable(double *x, const int n, const double rho, double (*oracle)(double),
        double l, double u);
void prox_l1(double *x, const int n, const double rho);
void prox_norm2(double *x, const double rho, const int n);
void prox_sum_square(double *x, const int n, const double rho);
void prox_log_barrier(double *x, const int n, const double rho);
void prox_huber(double *x, const int n, const double rho);
void prox_linear(double *x, const int n, const double rho, const double *c);
void prox_lp(double *x, const int n, const double rho, const double *c);

#endif
