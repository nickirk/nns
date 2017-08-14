#include <stdio.h>
#include <vector>
#include <random>
#include <math.h>
#include <algorithm>
#include <stdio.h>
#include <iostream>
#include "Basis.hpp"
#include "Determinant.hpp"
#include "detType.hpp"


Basis::Basis(int numOrb_, int numEle_){
  numEle = numEle_;
  numOrb = numOrb_; 
  indexOfDet = 0;
  size = calcSize(numOrb_, numEle_);
  for (int i= 0; i < numOrb; ++i) {listOfOrbNum.push_back(i);}
  createBasisDet(0, numEle_);
}

detType Basis::getDetByIndex(int index){
  return basis[index];
}

int Basis::getIndexByDet(detType const & det_){
  int pos = std::find(basis.begin(), basis.end(), det_)-basis.begin();
  return indexBasis[pos];
}

int Basis::calcSize(int numOrb_, int numEle_){
  if (numEle_==0) return 1;
  return (numOrb_ * calcSize(numOrb_ -1, numEle_-1)) / numEle_;
}
int Basis::getSize() const {return size;}

void Basis::createBasisDet(int offset, int numEle_){
  if (numEle_ == 0 ){
    detType tmpDet(numOrb);
    for (int i = 0; i < combination.size(); ++i){
      create(tmpDet,combination[i]);
    }
    indexBasis.push_back(indexOfDet++);
    //std::cout << "index= " << indexOfDet << std::endl;
    //std::cout << "intCast= " << tmpDet.intCast() << std::endl;
    basis.push_back(tmpDet);
    return;
  }
  for (int i = offset; i <= listOfOrbNum.size() - numEle_; ++i){
    combination.push_back(listOfOrbNum[i]);
    createBasisDet(i+1, numEle_-1);
    combination.pop_back();
  }
}

