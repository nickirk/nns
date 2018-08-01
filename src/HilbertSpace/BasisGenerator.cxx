/*
 * FermionBasis.cxx
 *
 *  Created on: Jul 26, 2018
 *      Author: Liao, guther
 */

#include "BasisGenerator.hpp"
#include "../Hamiltonian/Hamiltonian.hpp"
#include "../utilities/Errors.hpp"
#include "../HilbertSpace/Determinant.hpp"
#include <vector>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <iostream>

namespace networkVMC {

BasisGenerator::BasisGenerator(SpinConfig const &sC_):sC(sC_),numOrb(sC_.numSpinOrbs()),
		listOfOrbNum(sC_.numSpinOrbs()),combination(0){
	// create a list of orbital indices (just the numbers from 0 to #orbs
	std::iota(listOfOrbNum.begin(),listOfOrbNum.end(),0);
}

//---------------------------------------------------------------------------------------------------//

Basis BasisGenerator::generateBasis(Hamiltonian const &H) const{
	std::vector<detType> baseVec;
	switch(H.type()){
	case(Heisenberg):
		baseVec = generateSpinBasis();
		break;
	default:
		baseVec = generateFermionBasis();
	}

	return Basis(baseVec);
}

//---------------------------------------------------------------------------------------------------//

Basis BasisGenerator::generateBasis() const{
	return Basis(generateFermionBasis());
}

//---------------------------------------------------------------------------------------------------//

std::vector<detType> BasisGenerator::generateFermionBasis() const{
  int numEle = sC(-1)+sC(1);
  std::vector<detType> tmp(0),basis(0);
  // create the basis of all determinants with numEle electrons
  createBasisDets(tmp,0,numEle);
  int numSpinUp{0};
  int numSpinDown{0};
  // now, pick those with the correct spin
  for (size_t s(0); s<tmp.size(); s++){
    numSpinDown=0;
    numSpinUp=0;
    for (int i(0);i<(numOrb/2); ++i){
    	numSpinUp+= (tmp[s][2*i])?1:0;
    	numSpinDown+= (tmp[s][2*i+1])?1:0;
    }
    if ((numSpinUp == sC(1)) & (numSpinDown == sC(-1))){
      basis.push_back(tmp[s]);
    }
  }
  return basis;
}

//---------------------------------------------------------------------------------------------------//

std::vector<detType> BasisGenerator::generateSpinBasis() const{
	// create all determinants with a fixed number of 1-entries, but no further restrictions
	// this can be used as a spin basis, or a basis of spinless fermions (they have exactly the same dets)
	int numSpinUps = sC(1);
	std::vector<detType> basis(0);
	// recursively create the basis dets
	createBasisDets(basis,0,numSpinUps);
	return basis;
}

//---------------------------------------------------------------------------------------------------//

void BasisGenerator::createBasisDets(std::vector<detType> &basis, int offset, int numEle) const{
  if (numEle == 0 ){
    detType tmpDet(listOfOrbNum.size(),0);
    for (size_t i = 0; i < combination.size(); ++i){
      create(tmpDet,combination[i]);
    }
    basis.push_back(tmpDet);
    return;
  }
  for (size_t i = offset; i <= listOfOrbNum.size() - numEle; ++i){
    combination.push_back(listOfOrbNum[i]);
    createBasisDets(basis,i+1, numEle-1);
    combination.pop_back();
  }
}

//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//

int calcBasisSize(int numOrb_, int numEle_){
  if (numEle_==0) return 1;
  return (numOrb_ * calcBasisSize(numOrb_ -1, numEle_-1)) / numEle_;
}


} /* namespace networkVMC */
