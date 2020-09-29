#include "AmericanDualRogers.hpp"
#include "Product.hpp"
#include "BSModel.hpp"

double AmericanDualRogers::MinimizeGradientDescent(double init, double strike, int nbSamples,int n, gsl_rng* rng) {
    double lambda = init;
    for (uint i=1; i<=n ; ++i) {
        lambda = lambda - 1.0 / (i+1) * Gradient(lambda, strike, nbSamples, 0.0001, rng);
    }
    
    return lambda;
}

double AmericanDualRogers::Gradient(double lambda, double strike, int nbSamples, double epsilon, gsl_rng* rng) {
    double T = _product->_nbTimeSteps * _product->_delta;
    double sigma = gsl_vector_get(_model->_sigma, 0);
    double sum = 0;
    double var = 0;
    gsl_matrix* asset = gsl_matrix_alloc(_product->_nbTimeSteps +1, _model->_size);
    gsl_vector* martingale = gsl_vector_alloc(_product->_nbTimeSteps +1);
    for (int i=0; i<nbSamples; i++) {
        _model->Simul_Path(asset, T, _product->_nbTimeSteps, rng);
        GetMartingale(martingale, strike, asset);
        double max = (GetMax(martingale, asset, strike, _product->_nbTimeSteps, lambda+epsilon) - GetMax(martingale, asset, strike, _product->_nbTimeSteps, lambda-epsilon) ) / (2*epsilon);
        sum+= max;
    }
    sum = sum/nbSamples;
    return sum;
}

void AmericanDualRogers::Price(double &prix, double &ic, double lambda,double strike, int nbSamples, gsl_rng* rng) {
    double T = _product->_nbTimeSteps * _product->_delta;
    double sigma = gsl_vector_get(_model->_sigma, 0);
    double sum = 0;
    double var = 0;
    gsl_matrix* asset = gsl_matrix_alloc(_product->_nbTimeSteps +1, _model->_size);
    gsl_matrix* asset_minus = gsl_matrix_alloc(_product->_nbTimeSteps +1, _model->_size);
    gsl_vector* martingale = gsl_vector_alloc(_product->_nbTimeSteps +1);
    gsl_vector* martingale_minus = gsl_vector_alloc(_product->_nbTimeSteps +1);
    for (int i=0; i<nbSamples; i++) {
        _model->Simul_Path(asset, T, _product->_nbTimeSteps, rng, asset_minus);
        GetMartingale(martingale, strike, asset);
        GetMartingale(martingale_minus, strike, asset_minus);
        double max = GetMax(martingale, asset, strike, _product->_nbTimeSteps, lambda);
        double max_minus = GetMax(martingale_minus, asset_minus, strike, _product->_nbTimeSteps, lambda);
        double max_mean= (max+max_minus) / 2;
        sum+=max_mean;
        var += max_mean*max_mean;
    }
    var = var/nbSamples;
    prix = sum/nbSamples;
    double varianceEstimator = var - pow(prix, 2);
    ic = 1.96 * sqrt(varianceEstimator) / sqrt(nbSamples);
}

void AmericanDualRogers::GetMartingale(gsl_vector* martingale, double strike, gsl_matrix* path) {
    double T = _product->_nbTimeSteps * _product->_delta;
    bool isPayoffPositive = false;
    if (_product->Payoff(0, path) > 0) {
        isPayoffPositive=true;
    } 
    gsl_vector_set(martingale, 0, 0.0);
    for (int k=1; k<= _product->_nbTimeSteps; k++) {
        if (isPayoffPositive == true) {
            gsl_vector_set(martingale, k, gsl_vector_get(martingale, k-1)
                    +exp(-_model->_r *k*_product->_delta) * _product->TheoriticalPrice(k*_product->_delta, T, path)
                    -exp(-_model->_r*(k-1)*_product->_delta) * _product->TheoriticalPrice((k-1)*_product->_delta, T, path));
        } else {
            gsl_vector_set(martingale, k, gsl_vector_get(martingale,k-1));
        }
        if (isPayoffPositive==false && _product->Payoff(k, path) > 0) {
            isPayoffPositive = true;
        }
    }
}

double AmericanDualRogers::GetMax(gsl_vector* martingale, gsl_matrix* asset, double strike, int nbTimeSteps, double lambda) {
    double max = _product->Payoff(0, asset) - lambda*gsl_vector_get(martingale, 0);
    for (uint k=1; k<= nbTimeSteps ; k++) {
        double Z_t = exp(-_model->_r*k*_product->_delta)*_product->Payoff(k, asset) - lambda* gsl_vector_get(martingale, k);
        if ( Z_t > max )
            max = Z_t;
    }
    
    return max;
}