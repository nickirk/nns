#ifndef Hamiltonian_DEFINED
#define Hamiltonian_DEFINED

#include "Determinant.hpp"
#include "Basis.hpp"

// type of the determiants
class Hamiltonian{
  public:
    Hamiltonian(int size_, double diagDelta_, double offDiagMagntd_, double
		offDiagNonZeroRatio_){}
    int getSize() const;
    int getSparseSize() const;
    void sparseAcess(int pos, int &row, int &col, double &value);
    double & operator () (detType i, detType j);
    double & operator () (int i, int j);
    void multiplyMatrixVector(double *v, double *w);
  private:
    int size;
    double diagDelta;
    double offDiagMagntd;
    double OffDiagNonZeroRatio;
    std::vector<int> rowVec, colVec;
    std::vector<double> valueVec;
    void initHamiltonian();
};

#endif
