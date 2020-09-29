#ifndef AMERICANDUALML_HPP
#define AMERICANDUALML_HPP
#include "Product.hpp"
#include "Model.hpp"
#include "MonteCarlo.hpp"

class AmericanDualML {
public:
    Model* _model;
    Product* _product;
    MonteCarlo* _monteCarlo;
    
    AmericanDualML(Model* model, Product* product, MonteCarlo* monteCarlo): 
    _model(model), _product(product), _monteCarlo(monteCarlo) {
        tempVect = gsl_vector_alloc(_model->_size);
        path_i = gsl_vector_alloc(_model->_size);
    }
    ~AmericanDualML() {
        gsl_vector_free(tempVect);
        gsl_vector_free(path_i);
        gsl_vector_free(tau);
        gsl_matrix_free(past_i);
    }
    void getStrategy(gsl_vector* tau, gsl_matrix* asset, int timeStep, double delta, double theta);
    void getMartingale(gsl_matrix* path, gsl_matrix* asset, gsl_vector* martingale, int timeStep, long K, double delta, double theta, gsl_rng* rng);
    void MLPrice(double &prix, double &std, double delta, double theta, int timeStep, long n0, long k0, long n1, int L, int kappa);
    double MaxComputer(gsl_vector* martingale, gsl_matrix* path, double theta, double delta, int timeStep);
private:
    gsl_vector* tempVect;
    gsl_matrix* past_i;
    gsl_vector* path_i;
    gsl_vector* tau;
};

#endif /* AMERICANDUALML_HPP */

