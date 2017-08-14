#include <stdio.h>
#include <vector>
#include <random>
#include <math.h>
#include <algorithm>
#include "Basis.hpp"
#include "Determinant.hpp"


Basis::Basis(int numEle_, int numOrb_){
  numEle = numEle_;
  numOrb = numOrb_; 
  indexOfDet = 0;
  size = calcSize(numOrb_, numEle_);
  for (int i= 0; i < numOrb; ++i) {listOfOrbNum.push_back(i);}
  createBasisDet(0, numEle_);
}

Determinant Basis::getDetByIndex(int const & index_){
  int pos = std::find(indexBasis.begin(), 
		  indexBasis.end(), index_) - indexBasis.begin();
  return basis[pos];
}

int Basis::getIndexByDet(Determinant const & det_){
  int pos = std::find(basis.begin(), basis.end(), det_)-basis.begin();
  return indexBasis[pos];
}

int Basis::calcSize(int numOrb_, int numEle_){
  if (numEle_==0) return 1;
  return (numOrb_ * calcSize(numOrb_ -1, numEle_-1)) / numEle_;
}
int Basis::getSize(){return size;}

void Basis::createBasisDet(int offset, int numEle_){
  if (numEle_ == 0 ){
    Determinant tmpDet(numOrb);
    for (int i = 0; i < combination.size(); ++i){
      tmpDet.create(combination[i]); 
    }
    indexBasis.push_back(indexOfDet++);
    basis.push_back(tmpDet);
    return;
  }
  for (int i = offset; i <= listOfOrbNum.size() - numEle; ++i){
    combination.push_back(listOfOrbNum[i]);
    createBasisDet(i+1, numEle_-1);
    combination.pop_back();
  }
}

