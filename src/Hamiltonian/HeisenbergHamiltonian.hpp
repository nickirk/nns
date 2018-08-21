/*
 * HeisenbergHamiltonian.hpp
 *
 *  Created on: Jul 26, 2018
 *      Author: guther
 */

#ifndef SRC_HAMILTONIAN_HEISENBERGHAMILTONIAN_HPP_
#define SRC_HAMILTONIAN_HEISENBERGHAMILTONIAN_HPP_

#include "LatticeHamiltonian.hpp"

namespace networkVMC {

class HeisenbergHamiltonian: public LatticeHamiltonian {
  public:
	template<typename ...Args>
	HeisenbergHamiltonian(double J_, Args ...args):
		LatticeHamiltonian(true, args...),J(J_){};
	template<typename ...Args>
	HeisenbergHamiltonian(double J_, bool pbc, Args ...args):
		LatticeHamiltonian(pbc, args...),J(J_){};

	double operator()(detType const &a, detType const &b) const;

	HType type() const {return Heisenberg;}
	virtual ~HeisenbergHamiltonian(){};
private:
	double J;
};

//---------------------------------------------------------------------------------------------------//

// convert an occupancy flag to the ms value (0 is down, 1 is up)
inline int msVal(bool spin){
	if(spin) return 1;
	return -1;
}

} /* namespace networkVMC */

#endif /* SRC_HAMILTONIAN_HEISENBERGHAMILTONIAN_HPP_ */
