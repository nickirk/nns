/*
 * Hamiltonian.hpp
 *
 *  Created on: Jul 26, 2018
 *      Author: guther
 */

#ifndef SRC_HAMILTONIAN_HAMILTONIAN_HPP_
#define SRC_HAMILTONIAN_HAMILTONIAN_HPP_

#include "../utilities/TypeDefine.hpp"

namespace networkVMC {

// generic base class for Hamiltonians

class Hamiltonian {
public:
	Hamiltonian();
	// Essentially, Hamiltonians need to be able to produce a matrix element given two determinants
	virtual double operator()(detType const &a, detType const &b) const =0;
	// remarkably, what they do not have, is information on the size of the input determinants
	// in general, Hamiltonians can treat arbitrary a,b, but specifications might constrain this
	virtual ~Hamiltonian();

	// this generates all states coupled to source
	// note that - in contrast to excitation generation - this is really
	// a property of the Hamiltonian
	virtual std::vector<detType> getCoupledStates(detType const &source) const = 0;

	// and then, they have a type
	virtual HType type() const{return Constant;};
};

} /* namespace networkVMC */

#endif /* SRC_HAMILTONIAN_HAMILTONIAN_HPP_ */
