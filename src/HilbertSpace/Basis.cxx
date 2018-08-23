/*
 * Basis.cxx
 * Copyright (c) 2017
 * Author: Ke Liao and Kai Guther
 * All rights reserved
 */

#include "../utilities/Errors.hpp"
#include "../Hamiltonian/Hamiltonian.hpp"
#include "Basis.hpp"
#include <algorithm>
#include "../utilities/RNGWrapper.hpp"

#include "BasisGenerator.hpp"

namespace networkVMC{

// use the BasisGenerator to set up the list of determinants
// either without the Hamiltonian (deprecated)
Basis::Basis(SpinConfig const &sC){
	BasisGenerator bGen(sC);
	*this = bGen.generateBasis();
}

// or while supplying the Hamiltonian
Basis::Basis(SpinConfig const &sC, Hamiltonian const &H){
	BasisGenerator bGen(sC);
	*this = bGen.generateBasis(H);
}

detType Basis::getDetByIndex(int index) const{
  // do some bound checking
  if(static_cast<size_t>(index)>=basis.size() || index < 0) throw OutOfRangeError(index);
  return basis[index];
}

int Basis::getIndexByDet(detType const & det_) const{
  // Look for the determinant in the FermionBasis
  auto pos = std::find(basis.begin(), basis.end(), det_);
  // and get it's index
  auto dist = std::distance(basis.begin(),pos);
  // If we did not find the determinant, something went wrong
  if(static_cast<size_t>(dist)==basis.size()) throw InvalidDeterminantError(det_);
  return dist;
}

//---------------------------------------------------------------------------//

// This creates some random determinant from a given Basis
detType getRandomDeterminant(Basis const &fullBasis){
  // Essentially, we just get a random number and then return the corresponding det
  RNGWrapper rng; // Wrapper class for Mersenne-Twister RNG
  auto randomDet = fullBasis.getDetByIndex(rng());
  return randomDet;
}

}

