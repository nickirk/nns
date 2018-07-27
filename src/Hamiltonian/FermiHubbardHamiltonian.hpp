/*
 * FermiHubbardHamiltonian.hpp
 *
 *  Created on: Jun 20, 2018
 *      Author: guther
 */

#ifndef SRC_HAMILTONIAN_FERMIHUBBARDHAMILTONIAN_HPP_
#define SRC_HAMILTONIAN_FERMIHUBBARDHAMILTONIAN_HPP_

#include "LatticeHamiltonian.hpp"

namespace networkVMC {

class FermiHubbardHamiltonian: public LatticeHamiltonian {
public:

	// generate an n-d lattice
	template<typename ...Args>
	FermiHubbardHamiltonian(double U_, double t_, Args ...args):
	LatticeHamiltonian(args...),U(U_),t(t_){};

	// get the hubbard matrix element
	double operator()(detType const &a, detType const &b) const;
	virtual ~FermiHubbardHamiltonian(){};

	// hubbard's type
	virtual HType type() const {return Hubbard;}
private:
	double U, t;
};

} /* namespace networkVMC */

#endif /* SRC_HAMILTONIAN_FERMIHUBBARDHAMILTONIAN_HPP_ */
