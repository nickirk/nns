/*
 * testRBMAIH.cxx
 *
 *  Created on: Aug 24, 2018
 *      Author: guther
 */

#include "testComponents.hpp"
#include "../src/NNWLib.hpp"

using namespace networkVMC;

int main(){
	// read in the Hamiltonian
	auto modelHam = readAbInitioHamiltonian("FCIDUMP",true);
	auto numOrbs = modelHam.getNumOrbs();
	int nUp = 3;
	int nDown = 3;
	SpinConfig sC(nUp,nDown,numOrbs);

	testRBMMetropolis(sC,modelHam);
}


