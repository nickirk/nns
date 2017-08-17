/*
 * Hamiltonian.hpp
 * Copyright (c) 2017
 * Author: Ke Liao and Kai Guther
 * All rights reserved
 */

#ifndef Hamiltonian_DEFINED
#define Hamiltonian_DEFINED
#include <vector>
#include "Determinant.hpp"
#include "Basis.hpp"


// type of the determiants
class Hamiltonian{
  public:
    Hamiltonian(double diagDelta_, double offDiagMagntd_, double
		offDiagNonZeroRatio_, Basis const & basis_);
    int getSize() const;
    int getSparseSize() const;
    void sparseAccess(int pos, int &row, int &col, double &value) const;
    double  operator () (detType const &i, detType const &j) const;
    double  operator () (int const i, int const j) const;
    void multiplyMatrixVector(double *v, double *w);
    //What are these two functions for?
    int lowerPos(int i) const;
    int upperPos(int i) const;
  private:
    int size;
    double diagDelta;
    double offDiagMagntd;

    double offDiagNonZeroRatio;
    std::vector<int> rowVec, colVec;
    std::vector<double> valueVec;
    void initHamiltonian();
    void sortH();
    double lookupValue(int const &i, int const &j);
    Basis const& basis;
};

#endif
