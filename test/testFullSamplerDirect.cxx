/*
 * testFullSamplerDirect.cxx
 *
 *  Created on: Jul 27, 2018
 *      Author: guther
 */

#include "defaultSystem.hpp"
#include "testComponents.hpp"

int main(){
	int numSites = 4;
	auto sC = generateDefaultSpinConfig(numSites);
	auto modelHam = generateDefaultHubbard(numSites);
	testDeterministicFullSampling(sC,modelHam);
}
