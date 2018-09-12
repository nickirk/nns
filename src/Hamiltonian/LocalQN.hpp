/*
 * LocalQN.hpp
 *
 *  Created on: Jul 30, 2018
 *      Author: guther
 */

#ifndef SRC_HAMILTONIAN_LOCALQN_HPP_
#define SRC_HAMILTONIAN_LOCALQN_HPP_

#include "../utilities/TypeDefine.hpp"

namespace networkVMC {

/**
 * \class LocalQN
 * Converts site indices to orbital indices using a local basis
 * (typically this is spin-up/down)
 */
class LocalQN {
  public:
	// concrete orbital and range of possible numbers (range = #states per site)
	/**
	 * \param[in] orbital global index of a single-particle state with the orbital to consider
	 * \param[in] range number of orbitals per site
	 */
	LocalQN(int orbital, int range_):val(orbital%range_),range(range_){};
	virtual ~LocalQN(){};

	/**
	 * Get the index of the orbital on the site with index siteIndex
	 * \param[in] siteIndex index of the site
	 * \return global index of orbital on this site
	 */
	int orbitalIndex(int siteIndex) const{
		return range*siteIndex+val;
	}

	/**
	 * Get the index of the lattice site corresponding a global orbital index
	 * \param[in] orbIndex global index of the orbital on some site
	 * \return spatial index of that site
	 */
	int spatialIndex(int orbIndex) const{
		return orbIndex/range;
	}
private:
	/// Site-local index of the orbital considered
	int val;
	/// Number of orbitals per site
	int range;
};

/**
 * \param[in] ht type of Hamiltonian
 * \return Number of spin-orbitals per site for that Hamiltonian type
 * This is an auxiliary function, we hopefully do not need on the long term,
 * if we construct the basis using the local basis size instead of Hamiltonian type
 *
 */
inline int mapTypeToRange(HType const &ht){
	// returns the number of bits needed per site
	switch(ht){
	// currently, this is simple: Heisenberg has 2-d local Hilberspace
	// everything else 4-d
	case Heisenberg:
		return 1;
	default:
		return 2;
	}
}

} /* namespace networkVMC */

#endif /* SRC_HAMILTONIAN_LOCALQN_HPP_ */
