/*
 * Basis.hpp
 * Copyright (c) 2017
 * Author: Ke Liao and Kai Guther
 * All rights reserved
 */

#ifndef Basis_DEFINED
#define Basis_DEFINED

#include <vector>
#include "Determinant.hpp"
#include "../utilities/SpinConfig.hpp"

namespace networkVMC{

class Hamiltonian;

// Basis class for converting integers to basis states and indexing orbitals
// It is essentially a map from determinants to indices and vice versa
class Basis{
  public:
	// the hamiltonians type determines how the basis is set up
	explicit Basis(SpinConfig const &sC);
    Basis(SpinConfig const &sC, Hamiltonian const &H);
// total size of the many-body basis
    std::size_t size() const {return basis.size();}
// Return the determinant with index 'index'
    detType getDetByIndex(int index) const; // can throw an OutOfRangeError
// Return the index of the determinant 'det_'
    int getIndexByDet(detType const & det_) const; // can throw an InvalidDeterminantError

	// Return the alpha/beta spin distribution
	SpinConfig const& getSpinConfig() const {return spinConfig;};
    virtual ~Basis(){};
    // determinants need to be accessed often, no need to add some
    // overhead by using an extra class here, better use an alias
  private:
    // the conserved quantum numbers (aka spin configuration)
    SpinConfig spinConfig;
    // the vector of determinants
    std::vector<detType > basis;
};

}

#endif
