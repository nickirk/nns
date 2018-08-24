/*
 * testTrainTrial.cxx
 *
 *  Created on: Aug 24, 2018
 *      Author: guther
 */

#include "testComponents.hpp"
#include "../src/NNWLib.hpp"

using namespace networkVMC;

int main(){
	double J = -1.0;
	int Lx = 3;
	int Ly = 3;
	HeisenbergHamiltonian modelHam(J,Lx,Ly);

	int numSites = Lx*Ly;
	int nUp = numSites/2;
	int nDown = numSites - nUp;
	SpinConfig sC(nUp,nDown,numSites);

	testTrialMetropolis(sC, modelHam);
}



