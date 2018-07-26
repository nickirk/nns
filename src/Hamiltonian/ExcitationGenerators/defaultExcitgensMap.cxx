/*
 * defaultExcitgensMap.cxx
 *
 *  Created on: Jun 21, 2018
 *      Author: guther
 */

#include "defaultExcitgensMap.hpp"

#include "../TwoBodyHamiltonian.hpp"
#include "AllExcitationGenerators.hpp"

namespace networkVMC{

// return an excitation generator that is suited as a default
// for a given Hamiltonian
excitgenPtr getDefaultExcitgen(Hamiltonian const &H,
		detType const &HF){
	excitgenPtr defaultExcitgen{nullptr};
	switch(H.type()){
	case Hubbard:
		defaultExcitgen = excitgenPtr{new RSHubbardExcitgen{}};
		break;
	case AbInitio:
		defaultExcitgen = excitgenPtr{new WeightedExcitgen{H,HF}};
		break;
	case Constant:
	case Heisenberg:
	default:
		defaultExcitgen = excitgenPtr{new UniformExcitgen{HF}};
	}

	return defaultExcitgen;
}

}


