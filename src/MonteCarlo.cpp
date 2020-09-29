#include "MonteCarlo.hpp"

void MonteCarlo::Price(double& prix, double& ic, int nbSamples, gsl_rng* rng,int maturityIndex, bool isEuler) {
    int nbTimeSteps = _product->_nbTimeSteps;
    if (isEuler == true) {
        nbTimeSteps = 5*_product->_nbTimeSteps;
    }
    double sum = 0;
    double var = 0;
    gsl_matrix* path = gsl_matrix_alloc(nbTimeSteps+1, _model->_size);
    gsl_matrix* asset = gsl_matrix_alloc(_product->_nbTimeSteps+1, _model->_size); 
    for (int i=0; i<nbSamples; i++) {
        _model->Simul_Path(path, _product->_delta * _product->_nbTimeSteps, nbTimeSteps, rng);
        GetPathFromEulerScheme(asset, path, nbTimeSteps);
        double payoffNum= _product->Payoff(maturityIndex, asset) / _product->Numeraire(maturityIndex, asset);
        sum += payoffNum;
        var += payoffNum*payoffNum;
    }
    var = var/nbSamples;
    prix = sum/nbSamples;
    double varianceEstimator = (var - pow(sum, 2));
    ic = 1.96 * sqrt(varianceEstimator) / sqrt(nbSamples);
}

void MonteCarlo::Price(double t, double& prix, double& ic, int nbSamples, gsl_matrix* past, gsl_rng* rng, int maturityIndex, bool isEuler) {
    int nbTimeSteps = _product->_nbTimeSteps;
    if (isEuler == true) {
        nbTimeSteps = 5*_product->_nbTimeSteps;
    }
    double sum = 0;
    double var = 0;
    gsl_matrix* path = gsl_matrix_alloc(nbTimeSteps+1, _model->_size);
    gsl_matrix* asset = gsl_matrix_alloc(_product->_nbTimeSteps+1, _model->_size); 
    for (int i=0; i<nbSamples; i++) {
        _model->Simul_Path(path, t, _product->_delta * _product->_nbTimeSteps, nbTimeSteps, past, rng);
        GetPathFromEulerScheme(asset, path, nbTimeSteps);
        //cout << gsl_matrix_get(asset, 0, 0) << endl;
        double payoffNum= _product->Payoff(maturityIndex, asset) / _product->Numeraire(maturityIndex, asset);
        sum += payoffNum;
        var += payoffNum*payoffNum;
    }
    var = var/nbSamples;
    prix = sum/nbSamples;
    double varianceEstimator = (var - pow(sum, 2));
    ic = 1.96 * sqrt(varianceEstimator) / sqrt(nbSamples);
}

void MonteCarlo::GetPathFromEulerScheme(gsl_matrix* path, gsl_matrix* eulerPath, int nbTimeSteps) {
    for(uint i=0; i<= _model->_size; i++) {
        gsl_matrix_get_row(tempVect, eulerPath, i*nbTimeSteps / _model->_size);
        gsl_matrix_set_row(path, i, tempVect);
    }
}