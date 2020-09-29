#ifndef LIBORMODEL_HPP
#define LIBORMODEL_HPP

class LiborModel : public Model {
public:
    double _g_infinity;
    double _a;
    double _b;
    double _c;
    double _phi;
    
    LiborModel(int size, gsl_vector* spots, double g_infinity, double a, double b, double c, double phi);
    // Libor by Log euler scheme
    void Simul_Path(gsl_matrix* path,double T, int nbTimeSteps, gsl_rng* rng, gsl_matrix* path_minus=nullptr);
    void Simul_Path(gsl_matrix* path, double t, double T, int nbTimeSteps, gsl_matrix* past, gsl_rng* rng, gsl_matrix* path_minus=nullptr);
    void GetPathFromEulerScheme(gsl_matrix* path, gsl_matrix* eulerPath, int nbTimeSteps);
    ~LiborModel() {
        gsl_vector_free(tempVect);
    }
    int GetCoarseIndexFromt(double t, double T);
    double g(double t);
private:
    gsl_vector* tempVect;
};

#endif /* LIBORMODEL_HPP */

