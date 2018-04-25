/*
 * Basis.hpp
 * Copyright (c) 2017
 * Author: Ke Liao and Kai Guther
 * All rights reserved
 */

#ifndef Basis_DEFINED
#define Basis_DEFINED

#include "Determinant.hpp"
#include "../utilities/SpinConfig.hpp"
//
// type of the determiants

// Basis class for converting integers to basis states and indexing orbitals
// It is essentially a map from determinants to indices and vice versa
class Basis{
  public:
    Basis(SpinConfig const &spinConfig_);
// total size of the many-body basis
    int getSize() const;
// Return the determinant with index 'index'
    detType getDetByIndex(int index) const;
// Return the index of the determinant 'det_'
    int getIndexByDet(detType const & det_) const;
// Return the alpha/beta spin distribution
    SpinConfig getSpinConfig() const {return spinConfig;};
  private:
    int numEle;
    SpinConfig spinConfig;
    int numOrb;
    int size;
    int indexOfDet;
    std::vector<int> listOfOrbNum;
    std::vector<int> combination;
// internal methods for handling the map
    int calcSize(int numOrb_, int numEle_);
    void createBasisDet(int offset, int numEle_);
    // determinants need to be accessed often, no need to add some
    // overhead by using an extra class here, better use an alias
    std::vector<detType > basis;
    std::vector<int> indexBasis;
};

#endif
