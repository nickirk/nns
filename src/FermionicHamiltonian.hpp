/*
 * FermionicHamiltonian.h
 *
 *  Created on: Nov 15, 2017
 *      Author: guther
 */

#ifndef SRC_FERMIONICHAMILTONIAN_HPP_
#define SRC_FERMIONICHAMILTONIAN_HPP_

#include "Hamiltonian.hpp"

// Implementation for hamiltonians with fermionic commutation relations
class FermionicHamiltonian: public Hamiltonian {
public:
	FermionicHamiltonian(int dimension):Hamiltonian(dimension){};
	virtual ~FermionicHamiltonian();
// The fermi sign is the main difference to the plain Hamiltonian
	int getFermiSign(detType const &alpha, int annihilatorIndex, int creatorIndex) const;
// These generate either a single random state coupled to source (with probability p)
    detType getRandomCoupledState(detType const &source, double &p) const;
// or all states coupled to source
    std::vector<detType> getCoupledStates(detType const &source) const;
};

FermionicHamiltonian generateFermiHubbard(int dim, double U, double t);

#endif /* SRC_FERMIONICHAMILTONIAN_HPP_ */
