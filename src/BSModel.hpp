#ifndef BSMODEL_HPP
#define BSMODEL_HPP
using namespace std;
#include <vector>
#include <gsl/gsl_linalg.h>
#include <random>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "gsl/gsl_cdf.h" 
#include "Model.hpp"

class BSModel : public Model {
public:
    double _r; //interest rate
    gsl_vector* _sigma; //volatilities
    
    BSModel(int size, double r, gsl_vector* sigma, gsl_vector* S0,
            gsl_matrix* correlations);
    
    void Simul_Path(gsl_matrix* path, double T, int nbTimeSteps, gsl_rng* rng, gsl_matrix* path_minus = nullptr);
    void Simul_Path(gsl_matrix* path, double t, double T, int nbTimeSteps, gsl_matrix* past, gsl_rng* rng, gsl_matrix* path_minus=nullptr);
    ~BSModel() {}
    std::vector<double> ConvertGslVtoVector(gsl_vector* vector, int size);
    void GetLowerMatrix(gsl_matrix* matrix, int nbLines, int nbColumns);
};


#endif /* BSMODEL_HPP */

