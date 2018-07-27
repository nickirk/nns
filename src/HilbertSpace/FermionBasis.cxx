/*
 * FermionBasis.cxx
 *
 *  Created on: Jul 26, 2018
 *      Author: Liao, guther
 */

#include "FermionBasis.hpp"
#include "../utilities/Errors.hpp"
#include "../HilbertSpace/Determinant.hpp"
#include <vector>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <iostream>

namespace networkVMC {

std::vector<detType> generateFermionBasis(SpinConfig const &spinConfig){
  int numEle = spinConfig(-1)+spinConfig(1);
  int numOrb = spinConfig.numSpinOrbs();
  int indexOfDet = 0;
  std::vector<int> listOfOrbNum(numOrb), combination(0);
  std::iota(listOfOrbNum.begin(),listOfOrbNum.end(),0);
  std::vector<detType> tmp(0),basis(0);
  createFermionBasisDet(tmp,combination,0,numEle,listOfOrbNum);
  int numSpinUp{0};
  int numSpinDown{0};
  for (size_t s(0); s<tmp.size(); s++){
    numSpinDown=0;
    numSpinUp=0;
    for (int i(0);i<(numOrb/2); ++i){
    	numSpinUp+= (tmp[s][2*i])?1:0;
    	numSpinDown+= (tmp[s][2*i+1])?1:0;
    }
    if ((numSpinUp == spinConfig(1)) & (numSpinDown == spinConfig(-1))){
      basis.push_back(tmp[s]);
    }
  }
  return basis;
}

int calcBasisSize(int numOrb_, int numEle_){
  if (numEle_==0) return 1;
  return (numOrb_ * calcBasisSize(numOrb_ -1, numEle_-1)) / numEle_;
}

void createFermionBasisDet(std::vector<detType> &basis, std::vector<int> &combination, int offset,
		int numEle, std::vector<int> const &listOfOrbNum){
  if (numEle == 0 ){
    detType tmpDet(listOfOrbNum.size(),0);
    for (size_t i = 0; i < combination.size(); ++i){
      create(tmpDet,combination[i]);
    }
    std::cout << "intCast= " << verbatimCast(tmpDet) << std::endl;
    basis.push_back(tmpDet);
    return;
  }
  for (size_t i = offset; i <= listOfOrbNum.size() - numEle; ++i){
    combination.push_back(listOfOrbNum[i]);
    createFermionBasisDet(basis,combination,i+1, numEle-1,listOfOrbNum);
    combination.pop_back();
  }
}

} /* namespace networkVMC */
