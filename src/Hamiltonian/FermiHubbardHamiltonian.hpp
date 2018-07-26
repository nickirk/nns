/*
 * FermiHubbardHamiltonian.hpp
 *
 *  Created on: Jun 20, 2018
 *      Author: guther
 */

#ifndef SRC_HAMILTONIAN_FERMIHUBBARDHAMILTONIAN_HPP_
#define SRC_HAMILTONIAN_FERMIHUBBARDHAMILTONIAN_HPP_

#include "FermionicHamiltonian.hpp"

namespace networkVMC {

class FermiHubbardHamiltonian: public FermionicHamiltonian {
public:
	FermiHubbardHamiltonian(int dimension):
		FermionicHamiltonian(dimension){
		// for fermi-hubbard, we never round matrix elements
		// this is controlled by the partExact flag, so disable it
		partExactFlag = true;
	};
	virtual ~FermiHubbardHamiltonian();
	// get all couplings to source
	virtual std::vector<detType> getCoupledStates(detType const &source) const;
	virtual HType type() const {return Hubbard;}
};

FermiHubbardHamiltonian generateFermiHubbard(int dim, double U, double t);

} /* namespace networkVMC */

#endif /* SRC_HAMILTONIAN_FERMIHUBBARDHAMILTONIAN_HPP_ */
