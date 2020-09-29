#include <bits/random.h>
#include <gsl/gsl_matrix.h>
#include <stdio.h>
#include "BSModel.hpp"

using namespace std;

BSModel::BSModel(int size, double r, gsl_vector* sigma, gsl_vector* S0,
            gsl_matrix* correlations): 
        Model(size, S0, correlations), _r(r), _sigma(sigma){}

void BSModel::Simul_Path(gsl_matrix* results, double T, int nbTimeSteps, gsl_rng* rng, gsl_matrix* path_minus) {
    gsl_matrix_set_row(results, 0, _spots);
    if ( path_minus != nullptr) 
        gsl_matrix_set_row(path_minus, 0, _spots);
    double scalarProduct = 0;

    for (int i=1; i<= nbTimeSteps; i++) {
        gsl_ran_multivariate_gaussian(rng, means, cov, gaussianVector);
        for (int d = 0; d< _size; d++) {
            double sigma_d = gsl_vector_get(_sigma, d);
            gsl_matrix_get_row(correlD, _correlations, d);
            gsl_blas_ddot(correlD, gaussianVector, &scalarProduct);
            double results_i_d = gsl_matrix_get(results, i-1, d) * exp((_r-pow(sigma_d,2)*0.5)*T/nbTimeSteps 
                    + sigma_d * sqrt(T/nbTimeSteps)*scalarProduct);
            gsl_matrix_set(results, i, d , results_i_d);
            if (path_minus != nullptr) {
                double path_minus_i_d = gsl_matrix_get(path_minus, i-1, d) * exp((_r-pow(sigma_d,2)*0.5)*T/nbTimeSteps 
                    - sigma_d * sqrt(T/nbTimeSteps)*scalarProduct);
                gsl_matrix_set(path_minus, i, d , path_minus_i_d);
            }
        }
    }
    
}

void BSModel::Simul_Path(gsl_matrix* results, double t, double T, int nbTimeSteps, gsl_matrix* past, gsl_rng* rng, gsl_matrix* path_minus) {
    int lastDatePast = (int) floor(t * (double) nbTimeSteps / T);
    gsl_vector* tempRow = gsl_vector_alloc(_size);
    for(int i = 0; i <= lastDatePast; i++) {
        gsl_matrix_get_row(tempRow, past, i);
        gsl_matrix_set_row(results, i, tempRow);
        if ( path_minus != nullptr) 
            gsl_matrix_set_row(path_minus, i, tempRow);
    }
    /// Simulation de la trajectoire
    double timeInterval = (lastDatePast + 1) * T / (double) nbTimeSteps - t;
    double sigma_d = 0;
    double spot = 0;
    double spot_minus = 0;
    double results_t_d = 0;
    double scalarProduct = 0;
    for (int i = lastDatePast + 1; i <= nbTimeSteps; i++) {
        gsl_ran_multivariate_gaussian(rng, means, cov, gaussianVector);
        for (int d = 0; d < _size; d++) {
            sigma_d = gsl_vector_get(_sigma, d);
            gsl_matrix_get_row(correlD, _correlations, d);
            gsl_blas_ddot(correlD, gaussianVector, &scalarProduct);
            spot = (i == lastDatePast + 1) ? gsl_matrix_get(past, past->size1 - 1, d) : gsl_matrix_get(results, i - 1, d);
            spot_minus = (i == lastDatePast + 1) ? gsl_matrix_get(past, past->size1 - 1, d) : gsl_matrix_get(path_minus, i - 1, d);
            results_t_d = spot * exp((_r - pow(sigma_d, 2) / 2) * timeInterval +
                                      sigma_d * sqrt(timeInterval) * scalarProduct);
            
            gsl_matrix_set(results, i, d, results_t_d);
            if (path_minus != nullptr) {
                double path_minus_i_d =  spot_minus * exp((_r - pow(sigma_d, 2) / 2) * timeInterval -
                                      sigma_d * sqrt(timeInterval) * scalarProduct);
                gsl_matrix_set(path_minus, i, d, path_minus_i_d);
            }
        }
        if (i == lastDatePast + 1) {
            timeInterval = T / nbTimeSteps;
        }
    }
    
}

std::vector<double> BSModel::ConvertGslVtoVector(gsl_vector* vector, int size) {
    std::vector<double> results;
    for(int i=0; i<size; i++) {
        results.push_back(gsl_vector_get(vector, i));
    }
    return results;
}

void BSModel::GetLowerMatrix(gsl_matrix* matrix, int nbLines, int nbColumns) {
    for(uint i=0; i<nbLines; i++) {
        for(uint j=i+1; j<nbColumns; j++) {
            gsl_matrix_set(matrix, i, j, 0);
        }
    }
}