#ifndef GEOMETRICMEANPUT_HPP
#define GEOMETRICMEANPUT_HPP
#include "BSModel.hpp"

class GeometricMeanPut : public Product {
public:
    double _strike;
    BSModel* _model;
    GeometricMeanPut(BSModel *model, double T, int nbTimeSteps, double strike): 
    Product(T/nbTimeSteps, nbTimeSteps, model->_size),_model(model), _strike(strike) {}
    
    double Payoff(int maturityIndex, gsl_matrix* path);
    double GetGeoMean(int index,gsl_matrix* path);
    double TheoriticalPrice(double t, double T, gsl_matrix* asset= nullptr);
    double GetGeoMeanSigma();
    double GetSumVol();
};

#endif /* GEOMETRICMEANPUT_HPP */

