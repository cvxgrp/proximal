/* ------------------------------------------------------------- */
/* proximal_matrix.c -- operators that need linear algebra       */
/*     by N. Parikh, E. Chu, S. Boyd                             */ 
/*                                                               */
/* For performance reasons, these methods do not do any error    */
/* checking (e.g., array length). By compiling with the -fopenmp */
/* flag, these will automatically parallelize across all         */
/* axailable cores. The argument order is typically the point    */
/* at which to evaluate, the parameter rho, and then any other   */
/* parameters that are part of the definition of f. Where        */
/* relexant, the functions should accept INFINITY arguments.     */
/* ------------------------------------------------------------- */

#include "proximal_matrix.h"

// ==================== PROJECTIONS ====================

// f = I(X >= 0)
void project_sdc(gsl_matrix *X) 
{
    int n = X->size1;
    gsl_vector *d = gsl_vector_alloc(n);
    gsl_matrix *V = gsl_matrix_alloc(n,n);

    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(n);
    gsl_eigen_symmv(X, d, V, w);
    gsl_eigen_symmv_free(w);

    int i;
    double d_i;
    gsl_matrix_set_zero(X);
    for (i = 0; i < n; i++) {
        d_i = gsl_vector_get(d, i);
        if (d_i > 0) {
            gsl_vector_view V_i = gsl_matrix_column(V, i);
            gsl_blas_dsyr(CblasLower, d_i, &V_i.vector, X);
        }
    }

    gsl_vector_free(d);
    gsl_matrix_free(V);
}

// f = I(||X||_2 <= 1)
void project_spectral_norm_ball(gsl_matrix *X) 
{
    gsl_matrix *V = gsl_matrix_alloc(X->size1, X->size2);
    gsl_vector *d = gsl_vector_alloc(X->size2);
    gsl_vector *tmp = gsl_vector_alloc(X->size2);
    gsl_linalg_SV_decomp(X, W, d, tmp);

    int i;
    double d_i;
    gsl_matrix_set_zero(X);
    for (i = 0; i < X->size2; i++) {
        d_i = fmax(1, gsl_vector_get(d, i));
        gsl_vector_view U_i = gsl_matrix_column(X, i);
        gsl_vector_view V_i = gsl_matrix_column(V, i);
        gsl_blas_dger(d_i, &U_i.vector, &V_i.vector, X);
    }

    gsl_vector_free(d);
    gsl_matrix_free(V);
    gsl_vector_free(tmp);
}

// ==================== OTHER ====================

// f = -log(det(X))
void prox_neg_log_det(gsl_matrix *X, const double rho) 
{
    int n = X->size1;
    gsl_vector *d = gsl_vector_alloc(n);
    gsl_matrix *V = gsl_matrix_alloc(n,n);

    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(n);
    gsl_eigen_symmv(X, d, V, w);
    gsl_eigen_symmv_free(w);

    int i;
    double d_i;
    gsl_matrix_set_zero(X);
    for (i = 0; i < n; i++) {
        d_i = gsl_vector_get(d, i);
        d_i = (d_i + sqrt(pow(d_i,2) + 4*rho))/(2*rho);
        gsl_vector_set(d, i, d_i);
        gsl_vector_view V_i = gsl_matrix_column(V, i);
        gsl_blas_dsyr(CblasLower, d_i, &V_i.vector, X);
    }

    gsl_vector_free(d);
    gsl_matrix_free(V);
}

// f = (1/2) x^T Ax + b^T x
void prox_quad(gsl_vector *x, const double rho, gsl_matrix *A, gsl_matrix *b) 
{
    gsl_matrix *I = gsl_matrix_alloc(A->size1);
    gsl_matrix_set_identity(I);
    gsl_matrix_scale(I, rho);
    gsl_matrix_add(I, A);

    gsl_vector_scale(x, rho);
    gsl_vector_scale(b, -1);
    gsl_vector_add(b, x);

    gsl_linalg_cholesky_decomp(I);
    gsl_linalg_cholesky_solve(I, b, x);

    gsl_matrix_free(I);
}
