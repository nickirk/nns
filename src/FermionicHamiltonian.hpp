/*
 * FermionicHamiltonian.h
 *
 *  Created on: Nov 15, 2017
 *      Author: guther
 */

#ifndef SRC_FERMIONICHAMILTONIAN_HPP_
#define SRC_FERMIONICHAMILTONIAN_HPP_

#include "Hamiltonian.hpp"

class FermionicHamiltonian: public Hamiltonian {
public:
	FermionicHamiltonian(int dimension):Hamiltonian(dimension){};
	virtual ~FermionicHamiltonian();
	int getFermiSign(detType const &alpha, int annihilatorIndex, int creatorIndex) const;
    detType getRandomCoupledState(detType const &source, double &p) const;
    std::vector<detType> getCoupledStates(detType const &source) const;
};

FermionicHamiltonian generateFermiHubbard(int dim, double U, double t);

#endif /* SRC_FERMIONICHAMILTONIAN_HPP_ */
