#ifndef Hamiltonian_DEFINED
#define Hamiltonian_DEFINED

#include "Determinant.hpp"
#include "Basis.hpp"

// type of the determiants
class Hamiltonian{
  public:
    Hamiltonian(int size_, double diagDelta_, double offDiagMagntd_, double
		offDiagNonZeroRatio_);
    int getSize() const;
    int getSparseSize() const;
    void sparseAccess(int pos, int &row, int &col, double &value);
    double & operator () (detType i, detType j);
    double & operator () (int i, int j){return (*this)(detCast(i),detCast(j));}
    void multiplyMatrixVector(double *v, double *w);
    int lowerPos(int i){if(!sorted) sortH(); return std::find(rowVec.begin(),rowVec.end(),i);}
    int upperPos(int i){if(!sorted) sortH(); return std::find_end(rowVec.begin(),rowVec.end(),i);}
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
};

#endif
