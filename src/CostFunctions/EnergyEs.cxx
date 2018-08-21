/*
 * EnergyEs.cxx
 *
 *  Created on: Jun 26, 2018
 *      Author: guther
 */

#include "EnergyEs.hpp"
#include "EnergyCF.hpp"
#include "EnergyEsPreFetched.hpp"
#include "EnergyEsMarkov.hpp"

namespace networkVMC {
template <typename F, typename coeffType>
EnergyEs<F, coeffType>::EnergyEs(Hamiltonian const &H_, int numCons_):H(H_),
		worker(new EnergyCF<F, coeffType>(H, numCons_)), numCons(numCons_){
}

template <typename F, typename coeffType>
EnergyEs<F, coeffType>::~EnergyEs() {
}

template <typename F, typename coeffType>
void EnergyEs<F, coeffType>::setUpCF(SamplerType const &sT){
	// create the energy estimator belonging to the sampler-type
	switch(sT){
	case Markov:
		// for markov-type samplers, use EnergyEsMarkov
		worker = DeepCpyUniquePtr<EnergyCFBaseClass<F, coeffType>>(new EnergyEsMarkov(H,numCons));
		break;
	case PreFetched:
	default:
		// by default, use the pre-fetched version (i.e. the 'normal')
		worker = DeepCpyUniquePtr<EnergyCFBaseClass<F, coeffType>>(new EnergyEsPreFetched(H,numCons));
	}
}

} /* namespace networkVMC */
