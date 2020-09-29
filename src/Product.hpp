#ifndef PRODUCT_HPP
#define PRODUCT_HPP

class Product {
public:
    int _nbTimeSteps;
    int _size;
    double _delta;
    
    Product();
    
    Product(double delta, int nbTimeSteps, int size): _delta(delta), _size(size), _nbTimeSteps(nbTimeSteps) {}
    
    virtual ~Product() {}
    
    virtual double Payoff(int maturityIndex, gsl_matrix* path)=0;
    virtual double Numeraire(int maturityIndex, gsl_matrix* path) { return 1; }
    virtual double TheoriticalPrice(double t, double T, gsl_matrix* asset= nullptr)=0;
};

#endif /* PRODUCT_HPP */

