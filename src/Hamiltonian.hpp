#ifndef Hamiltonian_DEFINED
#define Hamiltonian_DEFINED

#include <vector>
#include "Determinant.hpp"
#include "Basis.hpp"


// type of the determiants
class Hamiltonian{
  public:
    Hamiltonian(int size_, double diagDelta_, double offDiagMagntd_, double
		offDiagNonZeroRatio_, Basis const & basis_);
    int getSize() const;
    int getSparseSize() const;
    void sparseAccess(int pos, int &row, int &col, double &value);
    double & operator () (detType i, detType j);
    double & operator () (int i, int j){return (*this)(basis.getDetByIndex(i),basis.getDetByIndex(j));}
    void multiplyMatrixVector(double *v, double *w);
    int lowerPos(int i){if(!sorted) sortH(); return std::find(rowVec.begin(),rowVec.end(),i);}
    int upperPos(int i){if(!sorted) sortH(); return std::find(rowVec.begin(),rowVec.end(),i);}
  private:
    int size;
    double diagDelta;
    double offDiagMagntd;
    double offDiagNonZeroRatio;
    std::vector<int> rowVec, colVec;
    std::vector<double> valueVec;
    void initHamiltonian();
    void sortH();
    bool sorted;
    double lookupValue(int const &i, int const &j);
    Basis const& basis;
};

#endif
