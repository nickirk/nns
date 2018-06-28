/*
 * ConnectionGenerator.cxx
 *
 *  Created on: Jun 28, 2018
 *      Author: guther
 */

#include "ConnectionGenerator.hpp"
#include "../ExcitationGenerator.hpp"
#include "../defaultExcitgensMap.hpp"

namespace networkVMC {

std::vector<detType> sampleConnections(Hamiltonian const &H, detType const &base,
		int numConnections, std::vector<double> &pGen){
	// create numConnections Determinants
	std::vector<detType> output(numConnections);
	// reset the size of pGen, this way we guarantee that it returns with the right size
	pGen.resize(numConnections);
	// create a default excitation generator for H
	auto excitgen = getDefaultExcitgen(H,base);

	// generate numConnections determinants
	for(int i = 0; i < numConnections; ++i){
		// for now, we allow recurrent excitations, but we will unbias with pGen
		output[i] = excitgen->generateExcitation(base,pGen[i]);
	}

	return output;
}

} /* namespace networkVMC */
