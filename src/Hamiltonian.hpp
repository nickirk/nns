#ifndef Hamiltonian_DEFINED
#define Hamiltonian_DEFINED

#include <vector>
#include "Determinant.hpp"
#include "Basis.hpp"


// type of the determiants
class Hamiltonian{
  public:
    Hamiltonian(double diagDelta_, double offDiagMagntd_, double
		offDiagNonZeroRatio_, Basis const &basis_);
    int getSize() const;
    int getSparseSize() const;
    //void sparseAcess(int pos, int &row, int &col, double &value);
    double & operator () (Determinant const &i, Determinant const &j);
    double & operator () (int const &i, int const &j);
    void multiplyMatrixVector(double *v, double *w);
  private:
    int size;
    double diagDelta;
    double offDiagMagntd;
    double offDiagNonZeroRatio;
    std::vector<int> rowVec, colVec;
    std::vector<double> valueVec;
    void initHamiltonian();
    double lookupValue(int const &i, int const &j);
    Basis const & basis;
};

#endif
