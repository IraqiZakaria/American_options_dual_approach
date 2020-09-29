#include <gsl/gsl_matrix_double.h>

#include "LiborModel.hpp"

LiborModel::LiborModel(int size, gsl_vector* spots, double g_infinity, double a, 
        double b, double c, double phi): _a(a), Model(size, spots),
        _b(b), _c(c), _g_infinity(g_infinity), _phi(phi) {
    this->_correlations= gsl_matrix_alloc(size,size);
    for (uint i=0; i<size; i++) {
        for(uint j=0; j<size; j++) {
            gsl_matrix_set(_correlations, i, j, exp(-phi* fabs(i-j)));
        }
    }
    gsl_linalg_cholesky_decomp(_correlations);
    GetLowerMatrix(_correlations, size, size);
    tempVect = gsl_vector_alloc(_size);
}

void LiborModel::Simul_Path(gsl_matrix* results, double T, int nbTimeSteps, gsl_rng* rng, gsl_matrix* path_minus) {
    double delta = T/_size;
    gsl_matrix_set_row(results, 0, _spots);
    double scalarProduct = 0;
    for (uint k=1; k<= nbTimeSteps; k++) {
        double tk = k*T/nbTimeSteps;
        gsl_ran_multivariate_gaussian(rng, means, cov, gaussianVector);
        int coarseIndex = GetCoarseIndexFromt(tk, T);
        coarseIndex = coarseIndex > 0 ? coarseIndex-1: coarseIndex;
        for (uint i = coarseIndex; i< _size ; i++) {
            double T_i = (i+1) * T/_size;
            double sigma_i = _c * g(T_i-tk);
            double sum =0;
            for (uint j= coarseIndex ; j<=i ; j++) {
                double T_j = (j+1) * T/_size;
                double Ltk_j = gsl_matrix_get(results, k-1, j);
                sum+= (delta* Ltk_j * sigma_i*_c* g(T_j-tk)*exp(-_phi*fabs(i-j))) / (1 + delta * Ltk_j); 
            }
            sum -= pow(sigma_i, 2) * 0.5;
            gsl_matrix_get_row(correlD, _correlations, i);
            gsl_blas_ddot(correlD, gaussianVector, &scalarProduct);
            double Ltk_i = gsl_matrix_get(results, k-1, i)*exp(sum*T/nbTimeSteps + sigma_i * sqrt(T/nbTimeSteps) * scalarProduct);
            gsl_matrix_set(results,k,i,Ltk_i);
        }
    }
}

void LiborModel::Simul_Path(gsl_matrix* results, double t, double T, int nbTimeSteps, gsl_matrix* past, gsl_rng* rng, gsl_matrix* path_minus) {
    int lastDatePast = (int) floor(t * (double) nbTimeSteps / T);
    double delta = T/_size;
    gsl_vector* tempRow = gsl_vector_alloc(_size);
    for(int i = 0; i <= lastDatePast; i++) {
        gsl_matrix_get_row(tempRow, past, i);
        gsl_matrix_set_row(results, i, tempRow);
    }
     //pnl_mat_get_row(spots_t, past, past->m - 1);

    /// Simulation de la trajectoire
    double timeInterval = (lastDatePast + 1) * T / (double) nbTimeSteps - t;
    //double sigma_d = 0;
    double spot = 0;
    //double results_t_d = 0;
    double scalarProduct = 0;
    for (uint k=lastDatePast + 1; k<= nbTimeSteps; k++) {
        double tk = k*T/nbTimeSteps;
        gsl_ran_multivariate_gaussian(rng, means, cov, gaussianVector);
        int coarseIndex = GetCoarseIndexFromt(tk, T);
        coarseIndex = coarseIndex > 0 ? coarseIndex-1: coarseIndex;
        for (uint i = coarseIndex; i< _size ; i++) {
            double T_i = (i+1) * T/_size;
            double sigma_i = _c * g(T_i-tk);
            double sum =0;
            for (uint j= coarseIndex ; j<=i ; j++) {
                double T_j = (j+1) * T/_size;
                double Ltk_j = gsl_matrix_get(results, k-1, j);
                sum+= (delta* Ltk_j * sigma_i*_c* g(T_j-tk)*exp(-_phi*fabs(i-j))) / (1 + delta * Ltk_j); 
            }
            sum -= pow(sigma_i, 2) * 0.5;
            gsl_matrix_get_row(correlD, _correlations, i);
            gsl_blas_ddot(correlD, gaussianVector, &scalarProduct);
            spot = (k == lastDatePast + 1) ? gsl_matrix_get(past, lastDatePast, i) : gsl_matrix_get(results, k - 1, i);
            double Ltk_i = spot*exp(sum*timeInterval + sigma_i * sqrt(timeInterval) * scalarProduct);
            gsl_matrix_set(results,k,i,Ltk_i);
        }
        if (k == lastDatePast + 1) {
            timeInterval = T / nbTimeSteps;
        }
    }
}

int LiborModel::GetCoarseIndexFromt(double t, double T) {
    return (int) floor(t * (double) _size / T);
}

double LiborModel::g(double t) {
    return _g_infinity + (1 - _g_infinity+ _a * t) * exp(-_b * t);
}

void LiborModel::GetPathFromEulerScheme(gsl_matrix* path, gsl_matrix* eulerPath, int nbTimeSteps) {
    for(uint i=0; i<= _size; i++) {
        gsl_matrix_get_row(tempVect, eulerPath, i*nbTimeSteps / _size);
        gsl_matrix_set_row(path, i, tempVect);
    }
}