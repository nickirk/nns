/*
 * BosonicHamiltonian.hpp
 *
 *  Created on: Nov 15, 2017
 *      Author: guther
 */

#ifndef SRC_BOSONICHAMILTONIAN_HPP_
#define SRC_BOSONICHAMILTONIAN_HPP_

#include "Hamiltonian.hpp"

// Hamiltonian with bosonic commutation relations
class BosonicHamiltonian: public Hamiltonian {
public:
	BosonicHamiltonian(int dimension):Hamiltonian(dimension){};
	virtual ~BosonicHamiltonian();
// Bosons dont have a Fermi sign, so we return 1
	int getFermiSign(detType const &alpha, int annihilatorIndex, int creatorIndex) const{return 1;};
    //detType getRandomCoupledState(detType const &source, double &p){};
    //std::vector<detType> getCoupledStates(detType const &source){} const;
};

#endif /* SRC_BOSONICHAMILTONIAN_HPP_ */
