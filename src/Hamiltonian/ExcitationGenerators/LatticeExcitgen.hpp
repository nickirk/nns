/*
 * LatticeExcitgen.hpp
 *
 *  Created on: Jul 30, 2018
 *      Author: guther
 */

#ifndef SRC_HAMILTONIAN_EXCITATIONGENERATORS_LATTICEEXCITGEN_HPP_
#define SRC_HAMILTONIAN_EXCITATIONGENERATORS_LATTICEEXCITGEN_HPP_

#include "ExcitationGenerator.hpp"
#include "../LatticeHamiltonian.hpp"

namespace networkVMC {

class LatticeExcitgen: public ClonableExcitgen<LatticeExcitgen> {
public:
	LatticeExcitgen(LatticeHamiltonian const &HL_):HL(HL_){};
	virtual ~LatticeExcitgen();

	// create an excitation on a lattice, i.e. using the adjacency list of the Hamiltonian
	detType generateExcitation(detType const &source, double &pGen);
	// get the excitation prob
	double getExcitationProb(detType const &source, detType const &target);
private:
	// underlying Hamiltonian, defining the lattice
	LatticeHamiltonian const &HL;
};

class LocalQN;
std::vector<int> getOccAdjs(int i, detType const &source, LatticeHamiltonian const &HL, LocalQN const &lQN);

} /* namespace networkVMC */

#endif /* SRC_HAMILTONIAN_EXCITATIONGENERATORS_LATTICEEXCITGEN_HPP_ */
