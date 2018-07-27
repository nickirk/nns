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

namespace networkVMC{

// Basis class for converting integers to basis states and indexing orbitals
// It is essentially a map from determinants to indices and vice versa
class Basis{
  public:
    Basis():size(0){};
// total size of the many-body basis
    std::size_t getSize() const {return size;}
// Return the determinant with index 'index'
    virtual detType getDetByIndex(int index) const = 0; // can throw an OutOfRangeError
// Return the index of the determinant 'det_'
    virtual int getIndexByDet(detType const & det_) const = 0; // can throw an InvalidDeterminantError
    virtual ~Basis(){};
  protected:
    std::size_t size;
};

}

#endif
