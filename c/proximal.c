/* ------------------------------------------------------------- */
/* proximal.c -- not slow implementations of proximal operators  */
/*     by N. Parikh, E. Chu, S. Boyd                             */ 
/*                                                               */
/* For performance reasons, these methods do not do any error    */
/* checking (e.g., array length). By compiling with the -fopenmp */
/* flag, these will automatically parallelize across all         */
/* available cores. The argument order is typically the point    */
/* at which to evaluate, the parameter rho, and then any other   */
/* parameters that are part of the definition of f. Where        */
/* relevant, the functions should accept INFINITY arguments.     */
/*                                                               */
/* Operators that require linear algebra to evaluate are in a    */
/* separate file (and require the GNU Scientific Library).       */
/* ------------------------------------------------------------- */

#include "proximal.h"

// ==================== PROJECTIONS ====================

// f = I(l <= x <= u)
void project_box(double *x, const int n, const double *l, const double *u) 
{
    int i;
#pragma omp parallel for
    for (i = 0; i < n; i++) {
        x[i] = fmax(l[i], fmin(x[i], u[i]));
    }
}

// f = I(x >= 0)
void project_pos(double *x, const int n)
{
    int i;
#pragma omp parallel for
    for (i = 0; i < n; i++) {
        x[i] = fmax(0, x[i]);
    }
}

// f = I(||x_{1:n}||_2 <= x_0)
void project_soc(double *x, const int n)
{
    double nx = 0;
    int i;

#pragma omp parallel for
    for (i = 1; i < n; i++) {
        nx += x[i]*x[i];
    }
    nx = sqrt(nx);

    if (nx <= x[0]) {
        // do nothing
    } else if (nx <= -x[0]) {
#pragma omp parallel for
        for (i = 0; i < n; i++) x[i] = 0;
    } else {
        double alpha = 0.5*(1 + x[0]/nx);
        x[0] = alpha*nx;
#pragma omp parallel for
        for (i = 1; i < n; i++) x[i] = alpha*x[i];
    }
}

// ==================== SEPARABLE FUNCTION ====================

double prox_scalar(double v, double rho, double (*oracle)(double), 
        double l, double u, double x0) 
{
    int MAX_ITER = 1000;
    double tol = 1e-8;
    double g;
    double x = fmax(l, fmin(x0, u));

    int iter;
    for (iter = 0; iter < MAX_ITER && u-l > tol; iter++) {
        g = -1/x + rho*(x - v);

        if (g > 0) {
            l = fmax(l, x - g/rho);
            u = x;
        } else if (g < 0) {
            l = x;
            u = fmin(u, x - g/rho);
        }

        x = (l + u)/2;
    }

    return x;
}

void prox_separable(double *x, const int n, const double rho, double (*oracle)(double),
        double l, double u)
{
    int i;
#pragma omp parallel for
    for (i = 0; i < n; i++) {
        x[i] = prox_scalar(x[i], rho, oracle, l, u, 0);
    }
}

// ==================== NORMS ====================

// f = ||.||_1
void prox_l1(double *x, const int n, const double rho) 
{
    int i;
#pragma omp parallel for
    for (i = 0; i < n; i++) {
        x[i] = fmax(0, x[i] - 1.0/rho) - fmax(0, -x[i] - 1.0/rho);
    }
}

// f = ||.||_2
void prox_norm2(double *x, const double rho, const int n) 
{
    double norm_square = 0;
    int i;

#pragma omp parallel for
    for (i = 0; i < n; i++) {
        norm_square += x[i]*x[i];
    }

    double norm = sqrt(norm_square);

#pragma omp parallel for
    for (i = 0; i < n; i++) {
        if (norm >= 1./rho)   x[i] *= (1. - 1./(rho*norm));
        else                  x[i] = 0;
    }
}

// ==================== OTHER ====================

// f = (1/2)||.||_2^2
void prox_sum_square(double *x, const int n, const double rho) 
{
    int i;
#pragma omp parallel for
    for (i = 0; i < n; i++) {
        x[i] *= rho/(1+rho);
    }
}

// f = -sum(log(x))
void prox_log_barrier(double *x, const int n, const double rho)
{
    int i;
#pragma omp parallel for
    for (i = 0; i < n; i++) {
        x[i] = 0.5*(x[i] + sqrt(x[i]*x[i] + 4./rho));
    }
}

// f = huber = x^2 if |x|<=1, 2|x| - 1 otherwise
double subgrad_huber(double x) {
    if (fabs(x) <= 1) {
        return 2*x;
    } else {
        return 2 * (x > 0 ? x : -x);
    }
}

void prox_huber(double *x, const int n, const double rho)
{
    prox_separable(x, n, rho, &subgrad_huber, -INFINITY, INFINITY);
}

// f = c'*x
void prox_linear(double *x, const int n, const double rho, const double *c) 
{
    int i;
#pragma omp parallel for
    for (i = 0; i < n; i ++){
        x[i] -= c[i]/rho;
    }
}

// f = c'*x + I(x >= 0)
void prox_lp(double *x, const int n, const double rho, const double *c) 
{
    int i;
#pragma omp parallel for
    for (i = 0; i < n; i ++){
        x[i] = fmax(0, x[i] - c[i]/rho);
    }
}
