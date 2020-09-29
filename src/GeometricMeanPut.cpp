#include "GeometricMeanPut.hpp"

double GeometricMeanPut::Payoff(int maturityIndex, gsl_matrix* path) {
    double geoMean = GetGeoMean(maturityIndex, path);
    return (_strike > geoMean ) ? (_strike - geoMean) : 0;
}

double GeometricMeanPut::GetGeoMean(int index, gsl_matrix* path) {
    double geoMean=1;
    for (uint i=0; i < _model->_size; i++) {
        geoMean *= gsl_matrix_get(path, index, i);
    }
    
    return pow(geoMean, 1.0 / _model->_size);
}

double GeometricMeanPut::TheoriticalPrice(double t, double T, gsl_matrix* asset) {
    double sigma = GetGeoMeanSigma();
    double sumSigmaSquare = GetSumVol();
    int t_index = round(t/_delta);
    double St = GetGeoMean(t_index, asset);
    double q = sumSigmaSquare / (2*_model->_size) - pow(sigma, 2) / 2 ;
    if (fabs(T-t) < 0.0001) {
        return Payoff(_nbTimeSteps, asset);
    }
    double d1 =(1.0/(sigma*sqrt(T-t)))*log(St/(_strike*exp(-(_model->_r - q)*(T-t)))) + sigma*sqrt(T-t)/2;
    double d0=d1 - sigma*sqrt(T-t);
    return _strike*exp(-_model->_r*(T-t))*gsl_cdf_gaussian_P(-d0, 1.0) - St*exp(-q*(T-t))*gsl_cdf_gaussian_P(-d1,1.0);
}

double GeometricMeanPut::GetGeoMeanSigma() {
    double sigma =0;
    for (uint i=0; i< _model->_size ; i++) {
        double sigma_i = gsl_vector_get(_model->_sigma, i);
        sigma += pow(sigma_i, 2);
        for (uint j= i+1; j< _model->_size ; j++) {
            double sigma_j = gsl_vector_get(_model->_sigma, j);
            double corr_ij = gsl_matrix_get(_model->_correlations, j, i);
            sigma += 2.0*sigma_i*sigma_j*corr_ij;
        }
    }
    
    return sqrt(sigma) / _model->_size;
}

double GeometricMeanPut::GetSumVol() {
    double sigma = 0;
    for (uint i=0; i< _model->_size ; i++) {
        sigma += pow(gsl_vector_get(_model->_sigma, i), 2);
    }
    
    return sigma;
}