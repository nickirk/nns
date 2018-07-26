/*
 * BosonicHamiltonian.hpp
 *
 *  Created on: Nov 15, 2017
 *      Author: guther
 */

#ifndef SRC_HAMILTONIAN_BOSONICHAMILTONIAN_HPP_
#define SRC_HAMILTONIAN_BOSONICHAMILTONIAN_HPP_

#include "TwoBodyHamiltonian.hpp"

namespace networkVMC{

// Hamiltonian with bosonic commutation relations
class BosonicHamiltonian: public TwoBodyHamiltonian {
public:
	BosonicHamiltonian(int dimension):TwoBodyHamiltonian(dimension){};
	virtual ~BosonicHamiltonian();
// Bosons dont have a Fermi sign, so we return 1
	int getFermiSign(detType const &alpha, int annihilatorIndex, int creatorIndex) const{return 0;};
    //detType getRandomCoupledState(detType const &source, double &p){};
    std::vector<detType> getCoupledStates(detType const &source) const;
};

}

#endif /* SRC_HAMILTONIAN_BOSONICHAMILTONIAN_HPP_ */
