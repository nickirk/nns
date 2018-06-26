/*
 * EnergyEs.cxx
 *
 *  Created on: Jun 26, 2018
 *      Author: guther
 */

#include "EnergyEs.hpp"
#include "EnergyEsPreFetched.hpp"
#include "EnergyEsMarkov.hpp"

namespace networkVMC {

EnergyEs::EnergyEs(Hamiltonian const &H_):H(H_) {
}

EnergyEs::~EnergyEs() {
}

std::unique_ptr<CostFunction> EnergyEs::setUpCF(SamplerType const &sT) const{
	// create the energy estimator belonging to the sampler-type
	switch(sT){
	case Markov:
		// for markov-type samplers, use EnergyEsMarkov
		return std::unique_ptr<CostFunction>(new EnergyEsMarkov(H));
		break;
	case PreFetched:
	default:
		// by default, use the pre-fetched version (i.e. the 'normal')
		return std::unique_ptr<CostFunction>(new EnergyEsPreFetched(H));
	}
}

} /* namespace networkVMC */
