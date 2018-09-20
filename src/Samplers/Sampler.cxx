/*
 * Sampler.cxx
 * Copyright (c) 2017
 * Author: Ke Liao and Kai Guther
 * All rights reserved
 */
#include "Sampler.hpp"

#include <vector>
#include <random>
#include <algorithm>
#include <Eigen/Dense>

#include "../Hamiltonian/ExcitationGenerators/defaultExcitgensMap.hpp"
#include "../HilbertSpace/Basis.hpp"
#include "../Hamiltonian/ExcitationGenerators/ExcitationGenerator.hpp"

namespace networkVMC{

// explicit constructor
Sampler::Sampler(ExcitationGenerator const &eG_, int numDets_):
	  excitGen(eG_.clone()),
	  numDets(numDets_){};

//---------------------------------------------------------------------------//

// construct the ExcitationGenerator implicitly from the Hamiltonian
Sampler::Sampler(Hamiltonian const &H_,
		  detType const &HF_, int numDets_):
			  excitGen(getDefaultExcitgen(H_,HF_).release()),
			  numDets(numDets_){};

//---------------------------------------------------------------------------//

// Just take any connected determinant and let the Excitation generator decide what is connected
detType Sampler::getRandomConnection(detType const &startingPoint, double &p) const{
	return excitGen->generateExcitation(startingPoint, p);
}

//---------------------------------------------------------------------------//

double Sampler::getConnectionProb(detType const &source, detType const &target) const{
	return excitGen->getExcitationProb(source,target);
}

}
