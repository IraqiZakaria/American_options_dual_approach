#ifndef MONTECARLO_HPP
#define MONTECARLO_HPP
#include "Model.hpp"
#include "Product.hpp"

class MonteCarlo {
public:
    Model* _model;
    Product* _product;
    MonteCarlo(Model* model, Product* product): _model(model), _product(product) {
        tempVect = gsl_vector_alloc(model->_size);
    }
    ~MonteCarlo() {
        gsl_vector_free(matrix_row);
        gsl_vector_free(tempVect);
    }
    void Price(double &prix, double &ic, int nbSamples, gsl_rng* rng, int maturityIndex , bool isEuler = false);
    void Price(double t, double &prix, double &ic, int nbSamples, gsl_matrix* past,gsl_rng* rng, int maturityIndex, bool isEuler = false);
    void GetPathFromEulerScheme(gsl_matrix* path, gsl_matrix* eulerPath, int nbTimeSteps);
private:
    gsl_vector* matrix_row;
    gsl_vector* tempVect;
};

#endif /* MONTECARLO_HPP */

