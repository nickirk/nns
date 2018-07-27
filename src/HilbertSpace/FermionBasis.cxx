/*
 * FermionBasis.cxx
 *
 *  Created on: Jul 26, 2018
 *      Author: Liao, guther
 */

#include "FermionBasis.hpp"
#include <vector>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <iterator>

namespace networkVMC {

FermionBasis::FermionBasis(SpinConfig const &spinConfig_):
	spinConfig(spinConfig_)
{
  numEle = spinConfig(-1)+spinConfig(1);
  numOrb = spinConfig.numSpinOrbs();
  indexOfDet = 0;
  for (int i= 0; i < numOrb; ++i) {
	  listOfOrbNum.push_back(i);
  }
  createFermionBasisDet(0, numEle);
  std::vector<detType> FermionBasisNew{0};
  int numSpinUp{0};
  int numSpinDown{0};
  for (size_t s(0); s<basis.size(); s++){
    numSpinDown=0;
    numSpinUp=0;
    for (int i(0);i<(numOrb/2); ++i){
    	numSpinUp+= (basis[s][2*i])?1:0;
    	numSpinDown+= (basis[s][2*i+1])?1:0;
    }
    if ((numSpinUp == spinConfig(1)) & (numSpinDown == spinConfig(-1))){
      FermionBasisNew.push_back(basis[s]);
    }
  }
  size = FermionBasisNew.size();
  basis=FermionBasisNew;
}

detType FermionBasis::getDetByIndex(int index) const{
  // do some bound checking
  if(static_cast<size_t>(index)>=basis.size() || index < 0) throw OutOfRangeError(index);
  return basis[index];
}

int FermionBasis::getIndexByDet(detType const & det_) const{
  // Look for the determinant in the FermionBasis
  auto pos = std::find(basis.begin(), basis.end(), det_);
  // and get it's index
  auto dist = std::distance(basis.begin(),pos);
  // If we did not find the determinant, something went wrong
  if(static_cast<size_t>(dist)==basis.size()) throw InvalidDeterminantError();
  return dist;
}

int FermionBasis::calcSize(int numOrb_, int numEle_){
  if (numEle_==0) return 1;
  return (numOrb_ * calcSize(numOrb_ -1, numEle_-1)) / numEle_;
}

void FermionBasis::createFermionBasisDet(int offset, int numEle_){
  if (numEle_ == 0 ){
    detType tmpDet(numOrb,0);
    for (size_t i = 0; i < combination.size(); ++i){
      create(tmpDet,combination[i]);
    }
    //std::cout << "intCast= " << tmpDet.intCast() << std::endl;
    basis.push_back(tmpDet);
    return;
  }
  for (size_t i = offset; i <= listOfOrbNum.size() - numEle_; ++i){
    combination.push_back(listOfOrbNum[i]);
    createFermionBasisDet(i+1, numEle_-1);
    combination.pop_back();
  }
}

} /* namespace networkVMC */
