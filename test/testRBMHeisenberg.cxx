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
	Basis basis(sC,HH);
	int const numHidden = 20;
	// RBM parametrization (needs complex coeffs)
	//RBM network(sC.numSpinOrbs(),numHidden);
	DirectParametrization<VecCType> network(basis);
	// start with the AFM determinant
	detType afmDet(sC.numSpinOrbs(),false);
	for(int i = 0; i < afmDet.size(); i+=2){
		afmDet[i] = true;
	}
	printDet(afmDet);
	// sampler
	//FullSampler<VecCType> mySampler(HH,basis,network);
	MetropolisSampler<VecCType> mySampler(HH,afmDet,network);
	solveADAM(network, mySampler, HH);
	//testRBMMetropolis(sC,HH);
}


