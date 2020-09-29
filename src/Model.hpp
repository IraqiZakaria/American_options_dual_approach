#ifndef MODEL_HPP
#define MODEL_HPP

#include <vector>
#include <gsl/gsl_linalg.h>
#include <random>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "gsl/gsl_cdf.h" 

class Model {
public:
    int _size; //number of underlyings in the model
    gsl_vector* _spots; //Spots
    gsl_matrix* _correlations; //Correlation matrix
    Model(int size, gsl_vector* spots);
    Model(int size, gsl_vector* spots, gsl_matrix* correlations);
    void GetLowerMatrix(gsl_matrix* matrix, int nbLines, int nbColumns);
    virtual void Simul_Path(gsl_matrix* path, double T, int nbTimeSteps, gsl_rng* rng, gsl_matrix* path_minus= nullptr)=0;
    virtual void Simul_Path(gsl_matrix* path, double t, double T, int nbTimeSteps, gsl_matrix* past, gsl_rng* rng, gsl_matrix* path_minus= nullptr)=0;
    virtual ~Model() {
        gsl_vector_free(means);
        gsl_matrix_free(cov);
        gsl_vector_free(correlD);
        gsl_vector_free(gaussianVector);
    };
protected:
    gsl_vector* means;
    gsl_matrix* cov;
    gsl_vector* correlD;
    gsl_vector* gaussianVector;
};

#endif /* MODEL_HPP */

