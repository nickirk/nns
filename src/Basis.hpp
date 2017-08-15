/*
 * Basis.hpp
 * Copyright (c) 2017
 * Author: Ke Liao and Kai Guther
 * All rights reserved
 */

#ifndef Basis_DEFINED
#define Basis_DEFINED

#include "Determinant.hpp"
#include "detType.hpp"
//
// type of the determiants



class Basis{
  public:
    Basis(int numEle_, int numOrb_);
    int getSize() const;
    detType getDetByIndex(int index) const;
    int getIndexByDet(detType const & det_) const;
  private:
    int numEle;
    int numOrb;
    int size;
    int indexOfDet;
    std::vector<int> listOfOrbNum;
    std::vector<int> combination;
    int calcSize(int numOrb_, int numEle_);
    void createBasisDet(int offset, int numEle_);
    // determinants need to be accessed often, no need to add some
    // overhead by using an extra class here, better use an alias
    std::vector<detType > basis;
    std::vector<int> indexBasis;
};

#endif
