#include "Model.hpp"

Model::Model(int size, gsl_vector* spots, gsl_matrix* correlations): 
        _size(size), _spots(spots){
    _correlations = gsl_matrix_alloc(size,size);
    gsl_matrix_memcpy(_correlations, correlations);
    gsl_linalg_cholesky_decomp(_correlations);
    GetLowerMatrix(_correlations, size, size);
    
    correlD = gsl_vector_alloc(_size);
    gaussianVector = gsl_vector_alloc(_size);
    cov = gsl_matrix_alloc(_size, _size);
    means = gsl_vector_alloc(_size);
    gsl_vector_set_zero(means);
    gsl_matrix_set_identity(cov);
}

Model::Model(int size, gsl_vector* spots): _size(size), _spots(spots) {
    correlD = gsl_vector_alloc(_size);
    gaussianVector = gsl_vector_alloc(_size);
    cov = gsl_matrix_alloc(_size, _size);
    means = gsl_vector_alloc(_size);
    gsl_vector_set_zero(means);
    gsl_matrix_set_identity(cov);
}

void Model::GetLowerMatrix(gsl_matrix* matrix, int nbLines, int nbColumns) {
    for(uint i=0; i<nbLines; i++) {
        for(uint j=i+1; j<nbColumns; j++) {
            gsl_matrix_set(matrix, i, j, 0);
        }
    }
}