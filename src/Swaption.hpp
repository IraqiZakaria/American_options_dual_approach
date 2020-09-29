#ifndef SWAPTION_HPP
#define SWAPTION_HPP
#include "Product.hpp"
#include "LiborModel.hpp"

class Swaption :public Product {
public:
    LiborModel* _model;
    double _theta;
    Swaption(LiborModel* model, double theta = 0.1, double delta = 0.25): _model(model), _theta(theta), Product(delta, model->_size, model->_size) {
        matrix_row = gsl_vector_alloc(model->_size);
        tempVect = gsl_vector_alloc(model->_size);
    }
    struct f_params {
        double c;
        double T1;
        double T2;
    };
    static double f(double t, void *p);
    //Approximation formula using Rebonato's formula
    double TheoriticalPrice(double t, double T, gsl_matrix* asset= nullptr);
    double Bond(double t, int maturityIndex, double delta, gsl_vector* L);
    double Payoff(int maturityIndex, gsl_matrix* L);
    double Numeraire(int maturityIndex, gsl_matrix* L);
    static double g(double t);
    void GetPathFromEulerScheme(gsl_matrix* path, gsl_matrix* eulerPath, int nbTimeSteps);
    ~Swaption(){
        gsl_vector_free(matrix_row);
        gsl_vector_free(tempVect);
    }
private:
    gsl_vector* matrix_row;
    gsl_vector* tempVect;

};

#endif /* SWAPTION_HPP */

