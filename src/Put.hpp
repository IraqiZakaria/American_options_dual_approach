#ifndef PUT_HPP
#define PUT_HPP

class Put : public Product {
public:
    double _strike;
    BSModel* _model;
    Put(BSModel* model, double T, int nbTimeSteps, double strike): _model(model), Product(T/nbTimeSteps, nbTimeSteps, 1), _strike(strike) {}
    double Payoff(int maturityIndex, gsl_matrix* path) {
        double S_T = gsl_matrix_get(path, maturityIndex, 0);
        return (_strike > S_T)? (_strike-S_T) : 0;
    }
    
    double TheoriticalPrice(double t, double T, gsl_matrix* asset= nullptr);
    
};

#endif /* PUT_HPP */

