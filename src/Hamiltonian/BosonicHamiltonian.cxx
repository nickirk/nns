/*
 * BosonicHamiltonian.cxx
 *
 *  Created on: Nov 15, 2017
 *      Author: guther
 */

#include "BosonicHamiltonian.hpp"

namespace networkVMC{

BosonicHamiltonian::~BosonicHamiltonian() {
	// TODO Auto-generated destructor stub
}

std::vector<detType> BosonicHamiltonian::getCoupledStates(
		detType const &source) const{
	return std::vector<detType>(0);
}

}

