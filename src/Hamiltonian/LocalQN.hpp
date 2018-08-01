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

// converts site indices to orbital indices using a local basis
// (typically this is spin-up/down)
class LocalQN {
public:
	// concrete orbital and range of possible numbers (range = #bits per site)
	LocalQN(int orbital, int range_):val(orbital%range_),range(range_){};
	virtual ~LocalQN(){};

	int orbitalIndex(int siteIndex) const{
		return range*siteIndex+val;
	}
	int spatialIndex(int orbIndex) const{
		return orbIndex/range;
	}
private:
	int val;
	int range;
};

// this is an auxiliary function, we hopefully do not need on the long term,
// if we construct the basis using the local basis size instead of Hamiltonian type
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
