/*
 * testRBMHeisenberg.cxx
 *
 *  Created on: Jul 30, 2018
 *      Author: guther
 */

#include "../src/NNWLib.hpp"
#include "testComponents.hpp"

using namespace networkVMC;

int main(){
	int numSites = 10;
	double J = -1.0;
	SpinConfig sC(numSites/2,numSites/2,numSites);
	HeisenbergHamiltonian HH(J,false,numSites);
	testAdj(HH);
	testDeterministicFullSampling(sC,HH);
}


