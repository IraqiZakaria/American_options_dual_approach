#include "Put.hpp"
#include "BSModel.hpp"

double Put::TheoriticalPrice(double t, double T, gsl_matrix* asset) {
    
    int index_t = round(t/ _delta);
    double sigma = gsl_vector_get(_model->_sigma, 0);
    double St = gsl_matrix_get(asset, index_t ,0);
    if (fabs(T-t) < 0.0001) {
        return Payoff(index_t, asset);
    }
    double d1 =(1.0/(sigma*sqrt(T-t)))*log(St/(_strike*exp(-_model->_r*(T-t)))) + sigma*sqrt(T-t)/2;
    double d0=d1 - sigma*sqrt(T-t);
    return _strike*exp(-_model->_r*(T-t))*gsl_cdf_gaussian_P(-d0, 1.0) - St*gsl_cdf_gaussian_P(-d1,1.0);
}