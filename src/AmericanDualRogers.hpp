#ifndef COMPUTATION_HPP
#define COMPUTATION_HPP
#include "Model.hpp"
#include "Product.hpp"
class AmericanDualRogers {
public:
    BSModel* _model;
    Product* _product;
    AmericanDualRogers() {}
    ~AmericanDualRogers() {
        gsl_matrix_free(path);
    }
    AmericanDualRogers(BSModel* model, Product* product): _model(model), _product(product) {}
    void Price(double &prix, double &ic, double lambda, double strike, int nbSamples, gsl_rng* rng);
    void GetMartingale(gsl_vector* martingales, double strike, gsl_matrix* path);
    double GetMax(gsl_vector* martingale, gsl_matrix* asset, double strike, int nbTimeSteps, double lambda);
    double Gradient(double lambda, double strike, int nbSamples, double epsilon, gsl_rng* rng);
    double GaussianCDF(double x, double mu, double sigma);
    double MinimizeGradientDescent(double lambda, double strike, int nbSamples, int n, gsl_rng* rng);
private:
    gsl_matrix* path;
};

#endif /* COMPUTATION_HPP */

