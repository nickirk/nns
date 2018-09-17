/*
 * testRBMHub.cxx
 *
 *  Created on: Sep 17, 2018
 *      Author: guther
 */

#include "../src/NNWLib.hpp"
#include "testComponents.hpp"

using namespace networkVMC;

int main(){
	double U = 4.;
	double t = -1.;
	int L = 10;
	int nStates = 2*L;
	FermiHubbardHamiltonian H(U,t,L);

	int nUp = L/2;
	int nDown = L/2 + L%2;
	SpinConfig sC(nUp,nDown,nStates);

	testRBMMetropolis(sC,H);


}


