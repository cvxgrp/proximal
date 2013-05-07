#ifndef PROXIMAL_MATRIX_H_GUARD
#define PROXIMAL_MATRIX_H_GUARD

#include <math.h>
#include <float.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

void project_sdc(gsl_matrix *X);
void project_spectral_norm_ball(gsl_matrix *X);
void prox_neg_log_det(gsl_matrix *X, const double rho);
void prox_quad(gsl_vector *x, const double rho, gsl_matrix *A, gsl_matrix *b);

#endif
