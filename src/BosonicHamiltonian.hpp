/*
 * BosonicHamiltonian.hpp
 *
 *  Created on: Nov 15, 2017
 *      Author: guther
 */

#ifndef SRC_BOSONICHAMILTONIAN_HPP_
#define SRC_BOSONICHAMILTONIAN_HPP_

#include "Hamiltonian.hpp"

class BosonicHamiltonian: public Hamiltonian {
public:
	BosonicHamiltonian(int dimension):Hamiltonian(dimension){};
	virtual ~BosonicHamiltonian();
	int getFermiSign(detType const &alpha, int annihilatorIndex, int creatorIndex) const{return 1;}
};

#endif /* SRC_BOSONICHAMILTONIAN_HPP_ */
