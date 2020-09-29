
#include "AmericanDualML.hpp"
#include "Model.hpp"
#include "Swaption.hpp"

void AmericanDualML::getStrategy(gsl_vector* tau, gsl_matrix* asset, int timeStep, double delta, double theta) {
    int lastMaturity = _model->_size/timeStep;
    for(uint i=0; i<= lastMaturity; i++) {
        for (uint j= i ; j<=lastMaturity; j++) {
            double max = 0; 
            double Z_j = _product->Payoff(j, asset) / _product->Numeraire(j, asset);
            for (uint p=j ; p<=lastMaturity; p++) {
                auto swaptionPrice = _product->TheoriticalPrice(j*delta, p*delta, asset);
                if ( swaptionPrice > max)
                    max = swaptionPrice;
                if ( swaptionPrice > Z_j)
                    break;
            }
            if (max <= Z_j) {
                gsl_vector_set(tau, i, j);
                break;
            }
        }
    }
}

void AmericanDualML::getMartingale(gsl_matrix* path, gsl_matrix* asset, gsl_vector* martingale,int timeStep, long K, double delta, double theta, gsl_rng* rng) {
    int nbTimeSteps = path->size1;
    double ic;
    int lastMaturity = _model->_size/timeStep;
    gsl_matrix_get_row(path_i, path, 0);
    gsl_matrix_set_row(past_i, 0, path_i);
    gsl_vector_set(martingale, 0, 0);
    getStrategy(tau, asset, timeStep, delta, theta);
    double priceBefore =0;
    gsl_matrix_get_row(path_i, path, 0);
    gsl_matrix_set_row(past_i, 0, path_i);
    bool tau_i_i = false;
    for (uint i= 1; i <= lastMaturity; i++) {
        int tau_i = gsl_vector_get(tau, i);
        double price;
        if (tau_i_i == true ) {
            price = -priceBefore;
            tau_i_i = false;
        } else {
            double prix;
            _monteCarlo->Price((i-1)*timeStep*delta, prix, ic, K, past_i, rng, tau_i, true);
            price = - prix;//_product->MonteCarloPricer(past_i, (i-1)*timeStep*delta, tau_i, delta, theta, K, rng);
        }
        for (uint j=(i-1)*timeStep*nbTimeSteps/_model->_size; j<=i*nbTimeSteps/_model->_size ; j++ ) {
                gsl_matrix_get_row(path_i, path, j);
                gsl_matrix_set_row(past_i, j, path_i);
        }
        if (tau_i <= i) {
            price += _product->Payoff(tau_i, path)/ _product->Numeraire(tau_i, path);
        } else {
            _monteCarlo->Price(i*timeStep*delta, priceBefore, ic, K, past_i, rng, tau_i, true);
            price += priceBefore;
            tau_i_i = true;
        }
        double martingale_i = gsl_vector_get(martingale, i-1);
        price += martingale_i;
        gsl_vector_set(martingale ,i , price);
    }
}

void AmericanDualML::MLPrice(double &prix, double &std, double delta, double theta, int timeStep, long n0, long k0, long n1, int L, int kappa) {
    const gsl_rng_type * Tt;
    gsl_rng * rng;

    gsl_rng_env_setup();

    Tt = gsl_rng_default;
    rng = gsl_rng_alloc (Tt);
    std::random_device rd{};
    int lastMaturity = _model->_size/timeStep;
    auto martingale = gsl_vector_alloc(lastMaturity+1);
    double sum=0;
    double var = 0;
    gsl_matrix* path = gsl_matrix_alloc(5*_model->_size+1, _model->_size);
    gsl_matrix* asset = gsl_matrix_alloc(_model->_size+1, _model->_size); 
    tau = gsl_vector_alloc(lastMaturity+1);
    past_i = gsl_matrix_alloc(5*_model->_size+1, _model->_size);
    for(uint i=0; i<n0; i++) {
         _model->Simul_Path(path, delta * _model->_size, 5*_model->_size, rng);
        _monteCarlo->GetPathFromEulerScheme(asset, path, 5*_model->_size);
        //cout << gsl_matrix_get(asset, 0, 0) << endl;
        getMartingale(path,asset, martingale, timeStep, k0, delta, theta, rng);
        double max= MaxComputer(martingale, path, theta , delta, timeStep);
        sum += max;
        var += max*max;
    }
    
    sum = sum / n0;
    var = var / n0;
    for(int l=1; l <=L ; l++) {
        int n_l = n1*pow(kappa, 1-l);
        int k_lMinus1 = k0 * pow(kappa, l-1);
        int k_l = k0 * pow(kappa, l);
        double sum_l = 0;
        double var_l = 0;
        for(uint i=0; i< n_l ; i++) {
             _model->Simul_Path(path, delta * _model->_size, 5*_model->_size, rng);
            _monteCarlo->GetPathFromEulerScheme(asset, path, 5*_model->_size);
            //cout << gsl_matrix_get(asset, 0, 0) << endl;
             getMartingale(path, asset, martingale, timeStep, k_l, delta, theta, rng);
            double max = MaxComputer(martingale, path, theta , delta, timeStep);
            getMartingale(path, asset, martingale, timeStep, k_lMinus1, delta, theta, rng);
            max= -MaxComputer(martingale, path, theta , delta, timeStep);
            sum_l += max;
            var_l += max*max;
        }
        sum +=  sum_l/n_l;
        var += var_l/n_l;
    }
    prix= sum;
    std = sqrt(var);
}

double AmericanDualML::MaxComputer(gsl_vector* martingale, gsl_matrix* path, double theta, double delta, int timeStep) {
    double Z_j = _product->Payoff(0, path) / _product->Numeraire(0, path);
    double max = Z_j- gsl_vector_get(martingale, 0);
    for (uint j=1; j< martingale->size ; j++) {
        double Z_j = _product->Payoff(j*timeStep, path) / _product->Numeraire(j* timeStep, path);
        double M_j = gsl_vector_get(martingale, j); 
        if ( Z_j - M_j > max )
            max = Z_j - M_j;
    }
    
    return max;
}