#include <iostream>
#include "BSModel.hpp"
#include "BSModel.cpp"
#include "Model.hpp"
#include "Model.cpp"
#include "LiborModel.hpp"
#include "LiborModel.cpp"
using namespace std;
#include <chrono>
#include "AmericanDualRogers.hpp"
#include "AmericanDualRogers.cpp"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "gsl/gsl_cdf.h" 
#include "Swaption.hpp"
#include "Swaption.cpp"
#include "AmericanDualML.hpp"
#include "AmericanDualML.cpp"
#include "Product.hpp"
#include "MonteCarlo.hpp"
#include "MonteCarlo.cpp"
#include "Put.hpp"
#include "Put.cpp"
#include "GeometricMeanPut.hpp"
#include "GeometricMeanPut.cpp"

int main(int argc, char**argv) {
    // Rogers 
    // 1. Put Américain
    int size=1; //number of underlyings in the model
    double r=0.06; //interest rate
    gsl_vector* sigma = gsl_vector_alloc(1); // vol
    gsl_vector_set_all(sigma, 0.4);
    gsl_matrix* correlations = gsl_matrix_alloc(1,1);//Correlation matrix
    gsl_matrix_set_all(correlations, 1);
    gsl_vector* spots = gsl_vector_alloc(1);//Spots
    gsl_vector_set_all(spots, 80);
    BSModel model(size, r, sigma, spots, correlations); //déclaration du modèle
    Put put(&model, 0.5, 50, 100); //déclaration de l'option (put)
    const gsl_rng_type * T;
    gsl_rng * rng;

    gsl_rng_env_setup();

    T = gsl_rng_default;
    rng = gsl_rng_alloc (T); //Initialisation du générateur
    AmericanDualRogers computer(&model, &put);
    double strike = 100;
    double lambda0 = 1;
    int M = 5000; //Nombre de trajectoires utilisé
    int nbTimeSteps = 50; //Nombre de timesteps utilisé
    int n = 300; //Nombre d'étapes utilisé pour l'algo du gradient stochastique
    std::cout << "Le prix du Put américain avec l'approche duale de Rogers:" << std::endl;
    std::cout << "strike =" << strike << std::endl;
    std::cout <<"Nombre de trajectoires= " << M << std::endl;
    std::cout <<"Nombre de timesteps=" << nbTimeSteps << std::endl;
    std::cout << "Volatilité= " << gsl_vector_get(sigma, 0) << std::endl;
    std::cout << "Taux d'intérêt=" << r << std::endl;
    std::cout << "Spot= " << gsl_vector_get(spots, 0) <<  std::endl;
    double prix,ic;
    //double lambda = computer.MinimizeGradientDescent(lambda, strike , M,n, rng);
    double lambda = 1.03913;
    computer.Price(prix, ic, lambda , strike, M, rng); // Calcul du prix ainsi que d'intervalle de confiance
   
    std::cout<< "Prix= " << prix << ", Intervalle de confiance=[" << prix - ic << "," << prix+ic<<"]" <<std::endl;
    std::cout <<"Lambda optimale=" << lambda <<std::endl;
    int sizeGeo = 2;
    double rho = 0.5;
    gsl_vector* sigmaGeo = gsl_vector_alloc(2); // vol
    gsl_vector_set(sigmaGeo,0, 0.4);
    gsl_vector_set(sigmaGeo,1, 0.1);
    gsl_matrix* correlationsGeo = gsl_matrix_alloc(2,2);//Correlation matrix
    gsl_matrix_set_all(correlationsGeo, 1);
    gsl_matrix_set(correlationsGeo, 0,1 , rho);
    gsl_matrix_set(correlationsGeo, 1, 0, rho);
    gsl_vector* spotsGeo = gsl_vector_alloc(2);//Spots
    gsl_vector_set_all(spotsGeo, 80);
    gsl_vector_set(spotsGeo, 1, 90);
    BSModel modelGeo(sizeGeo, r, sigmaGeo, spotsGeo, correlationsGeo); //déclaration du modèle
    GeometricMeanPut geoPut(&modelGeo, 0.5, 50, 100); //déclaration de l'option (put sur moyenne géométrique)
    AmericanDualRogers geoComputer(&modelGeo, &geoPut);
    double lambdageo = 1.07907;//geoComputer.MinimizeGradientDescent(lambda, strike , M,n, rng);
    geoComputer.Price(prix, ic, lambdageo , strike, M, rng);
    std::cout << "Le prix du put sur une moyenne géométrique" << std::endl;
    std::cout << "strike =" << strike << std::endl;
    std::cout <<"Nombre de trajectoires= " << M << std::endl;
    std::cout <<"Nombre de timesteps=" << nbTimeSteps << std::endl;
    std::cout << "1ère volatilité= " << gsl_vector_get(sigmaGeo, 0) << std::endl;
    std::cout << "2ème volatilité= " << gsl_vector_get(sigmaGeo, 1) << std::endl;
    std::cout << "Taux d'intérêt=" << r << std::endl;
    std::cout << "1er Spot= " << gsl_vector_get(spotsGeo, 0) <<  std::endl;
    std::cout << "2eme Spot= " << gsl_vector_get(spotsGeo, 1) <<  std::endl;
    std::cout << "prix= " << prix << " Intervalle de confiance=[" << prix - ic << "," << prix+ic<<"]"<< "lambda= " << lambdageo << std::endl;
    
    
    
    gsl_vector* S0 = gsl_vector_alloc(40);//Spots
    gsl_vector_set_all(S0, 0.1);
    LiborModel modsel(40, S0, 0.5, 1.5, 3.5 , 0.2, 0.0413);
    double theta = 0.1;
    double k0 = 50;
    double n0= 6600;
    double n1 = 2230;
    double std;
    Swaption swapchn(&modsel, theta, 0.25);
    std::cout << "Multilevel pricing(ça risque de prendre beaucoup de temps)" << std::endl;
    std::cout << "strike =" << theta << std::endl;
    std::cout <<"n0= " << n0 << std::endl;
    std::cout <<"n1= " << n1 << std::endl;
    std::cout <<"k0= " << k0 << std::endl;
    std::cout << "vaeur initiale= " << gsl_vector_get(S0, 0) <<  std::endl;
    MonteCarlo monteCarlo(&modsel, &swapchn);
    AmericanDualML multilevel(&modsel, &swapchn, &monteCarlo);
    multilevel.MLPrice(prix, std, 0.25, theta, 4, n0, k0, n1, 3, 2);
    std::cout << "prix= " << prix << "écart-type= " << std << std::endl;
    return 0;
}
