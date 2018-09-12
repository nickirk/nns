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

/**
 * \class LatticeExcitgen
 * \brief Generic excitationGenerator for lattice Hamiltonians
 *
 * This uses the adjacency defined on the Lattice to generate excitations. It is slightly
 * slower than the dedicated RSHubbardExcitgen for 1-d Hubbard systems, but much more efficient than
 * the other ExcitationGenerators when applied to lattice classes
 */
class LatticeExcitgen: public ClonableExcitgen<LatticeExcitgen> {
public:
	/**
	 * \param[in] HL_ hamiltonian defining the connectivity
	 */
	LatticeExcitgen(LatticeHamiltonian const &HL_):HL(HL_){};
	virtual ~LatticeExcitgen();

	// create an excitation on a lattice, i.e. using the adjacency list of the Hamiltonian
	detType generateExcitation(detType const &source, double &pGen);
	// get the excitation prob
	double getExcitationProb(detType const &source, detType const &target);
private:
	/// underlying Hamiltonian, defining the lattice
	LatticeHamiltonian const &HL;
};

class LocalQN;
/**
 * \fn std::vector<int> getOccAdjs(int i, detType const &source, LatticeHamiltonian const &HL, LocalQN const &lQN);
 * \brief Get all sites adjacent to a given one which are occupied in a given basis vector
 * \param[in] i lattice site to consider
 * \param[in] source basis state to consider
 * \param[in] HL hamiltonian defining connectivity
 * \param[in] lQN local basis state to consider
 * \return list of the indices of sites that have lQN at i occupied in source and adjacent to i
 */
std::vector<int> getOccAdjs(int i, detType const &source, LatticeHamiltonian const &HL, LocalQN const &lQN);

} /* namespace networkVMC */

#endif /* SRC_HAMILTONIAN_EXCITATIONGENERATORS_LATTICEEXCITGEN_HPP_ */
