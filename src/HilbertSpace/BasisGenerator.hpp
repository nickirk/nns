/*
 * FermionBasis.hpp
 *
 *  Created on: Jul 26, 2018
 *      Author: guther
 */

#ifndef SRC_HILBERTSPACE_BASISGENERATOR_HPP_
#define SRC_HILBERTSPACE_BASISGENERATOR_HPP_

#include "../utilities/SpinConfig.hpp"
#include <vector>
#include "../utilities/TypeDefine.hpp"
#include "Basis.hpp"

namespace networkVMC {

class Hamiltonian;

// This one generates basis lists
class BasisGenerator{
public:
	BasisGenerator(SpinConfig const &sC_);
	// create the basis matching H
	Basis generateBasis(Hamiltonian const &H)  const;
	// if no H is passed, create a fermion basis
	Basis generateBasis() const;
private:
	// the spin configuration (number of up/down spins)
	SpinConfig sC;
	// the number of spin-orbitals to occupy
	int numOrb;
	// a list of the orbitals
	std::vector<int> listOfOrbNum;
	// and a temporary for storing determinants as orbital lists
	mutable std::vector<int>combination;

	// the decision on which basis to take is made by the Hamiltonian
	// so we cannot call these directly
	std::vector<detType> generateFermionBasis() const;
	std::vector<detType> generateSpinBasis() const;

	// this recursively builds up a list of all determinants with the given number of electrons
	void createBasisDets(std::vector<detType> &basis, int offset, int numEle) const;
};

} /* namespace networkVMC */

#endif /* SRC_HILBERTSPACE_BASISGENERATOR_HPP_ */
