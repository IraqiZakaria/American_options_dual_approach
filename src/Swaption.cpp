#include "Swaption.hpp"
#include <gsl/gsl_integration.h>

double Swaption::Bond(double t, int maturityIndex, double delta, gsl_vector* L) {
    double result=1;
    for(uint i=_model->GetCoarseIndexFromt(t, _model->_size*delta)+1; i<=maturityIndex; i++) {
        result*=1. / (1+ gsl_vector_get(L, i) * delta);
    }
    
    return result;
}

double Swaption::TheoriticalPrice(double t, double T, gsl_matrix* asset) {
    int index = t/_delta;
    if ( t == T ) {
        return Payoff(T/_delta, asset) / Numeraire(T/_delta, asset);
    }
    double bondSum = 0.;
    double swapRate = 0;
    double volSquare = 0;
    gsl_vector* weights = gsl_vector_alloc(_model->_size);
    gsl_vector* L = gsl_vector_alloc(_model->_size);
    gsl_matrix_get_row(L, asset, index);
    f_params params;
    params.c=_model->_c;

    gsl_function F;
   

    double result, error;
    size_t neval;
    F.function = &f;
    const double epsabs=1e-4;
    const double epsrel=1e-4;
    int kappa_T = _model->GetCoarseIndexFromt(T, _model->_size*_delta);
    for (uint i=kappa_T; i< _model->_size; i++) {
        bondSum += Bond(0, i, _delta, L);
    }
    for(uint i=kappa_T; i< _model->_size ; i++) {
        double weight_i = Bond(0, i, _delta, L) / bondSum;
        gsl_vector_set(weights, i, weight_i);
    }
    for (uint i=kappa_T; i<_model->_size; i++) {
        swapRate += gsl_vector_get(weights, i)* gsl_vector_get(L, i);
    }
    
    for (uint i=kappa_T; i<_model->_size ; i++) {
        double w_i = gsl_vector_get(weights, i);
        double L_i = gsl_vector_get(L, i);
        double T_i = (i+1) * _delta;
        params.T1 = T_i;
        for(uint j=kappa_T; j < _model->_size; j++) {
            double T_j = (j+1) * _delta;
            params.T2 = T_j;
            F.params = reinterpret_cast<void *>(&params);
            double w_j = gsl_vector_get(weights, j);
            double L_j = gsl_vector_get(L, j);
            double rho_i_j = exp(-_model->_phi*fabs(i-j));
            gsl_integration_qng (&F, 
                                  t,
                                  T,
                                  epsabs, 
                                  epsrel, 
                                  &result,
                                  &error, 
                                  &neval);
            volSquare += (w_i*w_j*L_i*L_j*rho_i_j)/(swapRate*swapRate) * result;
        }
    }
    double d1 =(1.0/(sqrt(volSquare)))*log(swapRate/_theta) + sqrt(volSquare)/2;
    double d0=d1 - sqrt(volSquare);
    return _delta*bondSum*(swapRate*gsl_cdf_gaussian_P(d1,1.0) -  _theta*gsl_cdf_gaussian_P(d0, 1.0));
}

double Swaption::f(double t, void* p){
    f_params &params= *reinterpret_cast<f_params *>(p);
    return params.c * params.c* g(params.T1 - t) * g(params.T2 - t);
}

double Swaption::g(double t) {
    return 0.5 + (1 - 0.5+ 1.5 * t) * exp(-3.5 * t);
}

double Swaption::Payoff(int maturityIndex, gsl_matrix* L) {
    double payoff = 0;
    gsl_matrix_get_row(matrix_row ,L, maturityIndex);
    for(uint i= maturityIndex; i < _model->_size; i++) {
        double L_j_i = gsl_matrix_get(L, maturityIndex, i);
        payoff += Bond(maturityIndex * _delta, i, _delta, matrix_row)* _delta*(L_j_i - _theta);
    }
    return payoff > 0 ? payoff: 0;
}

double Swaption::Numeraire(int maturityIndex, gsl_matrix* L) {
    double result = 1;
    for (uint i=1; i< maturityIndex ; i++) {
        double L_j_i = gsl_matrix_get(L, i, i);
        result*= (1 + _delta * L_j_i);
    }
    gsl_matrix_get_row(matrix_row ,L, 0);
    result /= Bond(0, 0, _delta, matrix_row);
    return result;
}