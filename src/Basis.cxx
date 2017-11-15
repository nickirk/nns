/*
 * Basis.cxx
 * Copyright (c) 2017
 * Author: Ke Liao and Kai Guther
 * All rights reserved
 */

#include <stdio.h>
#include <vector>
#include <random>
#include <math.h>
#include <algorithm>
#include <stdio.h>
#include <iostream>
#include "Basis.hpp"
#include "Determinant.hpp"


Basis::Basis(std::vector<int> const &spinConfig){
  numEle = spinConfig[0]+spinConfig[1];
  numOrb = spinConfig[2]; 
  indexOfDet = 0;
  for (int i= 0; i < numOrb; ++i) {listOfOrbNum.push_back(i);}
  createBasisDet(0, numEle);
  std::vector<detType> basisNew;
  int numSpinUp(0);
  int numSpinDown(0);
  for (size_t s(0); s<basis.size(); s++){
    numSpinDown=0;
    numSpinUp=0;
    for (int i(0);i<(numOrb/2); ++i){
	numSpinUp+= (basis[s][2*i])?1:0;
	numSpinDown+= (basis[s][2*i+1])?1:0;
    }
    if ((numSpinUp == spinConfig[0]) & (numSpinDown == spinConfig[1])) 
      basisNew.push_back(basis[s]);
  }
  size = basisNew.size();
  basis=basisNew;
}

int Basis::getSize() const {return size;}

detType Basis::getDetByIndex(int index) const{
  return basis[index];
}

int Basis::getIndexByDet(detType const & det_) const{
  int pos = std::find(basis.begin(), basis.end(), det_)-basis.begin();
  return indexBasis[pos];
}

int Basis::calcSize(int numOrb_, int numEle_){
  if (numEle_==0) return 1;
  return (numOrb_ * calcSize(numOrb_ -1, numEle_-1)) / numEle_;
}

void Basis::createBasisDet(int offset, int numEle_){
  if (numEle_ == 0 ){
    detType tmpDet(numOrb,0);
    for (size_t i = 0; i < combination.size(); ++i){
      create(tmpDet,combination[i]);
    }
    indexBasis.push_back(indexOfDet++);
    //std::cout << "index= " << indexOfDet << std::endl;
    //std::cout << "intCast= " << tmpDet.intCast() << std::endl;
    basis.push_back(tmpDet);
    return;
  }
  for (size_t i = offset; i <= listOfOrbNum.size() - numEle_; ++i){
    combination.push_back(listOfOrbNum[i]);
    createBasisDet(i+1, numEle_-1);
    combination.pop_back();
  }
}

