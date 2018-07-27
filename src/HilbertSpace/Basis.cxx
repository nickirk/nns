/*
 * Basis.cxx
 * Copyright (c) 2017
 * Author: Ke Liao and Kai Guther
 * All rights reserved
 */

#include "../utilities/Errors.hpp"
#include "../Hamiltonian/Hamiltonian.hpp"
#include "Basis.hpp"
#include "FermionBasis.hpp"

namespace networkVMC{

Basis::Basis(SpinConfig const &sC):spinConfig(sC){
	basis = generateFermionBasis(sC);
}

Basis::Basis(SpinConfig const &sC, Hamiltonian const &H):spinConfig(sC){
	switch(H.type()){
	case(Heisenberg):
			break;
	default:
		basis = generateFermionBasis(sC);
	}
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

}

