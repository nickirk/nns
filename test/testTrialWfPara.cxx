/*
 * testTrialWfPara.cxx
 *
 *  Created on: Aug 16, 2018
 *      Author: guther
 */

#include <iostream>
#include <assert.h>
#include "../src/NNWLib.hpp"
#include "defaultSystem.hpp"

using namespace networkVMC;

int main(){
	int numSites = 4;
	int nHidden = 8;
	auto twfBasis = generateDefaultBasis(numSites);
	// do a simple test of the TrialWfPara with
	// an RBM as a base para
	RBM network(numSites*2,nHidden);
	auto HF = twfBasis.getDetByIndex(0);

	// and a linear parametrization as a trial WF
	auto trial = DirectParametrization<VecCType>(twfBasis);

	TrialWfPara<VecCType> twfPara(network,trial);

	// first compare the coefficients of the base*trial with that of the TrialWfPara
	assert(std::abs(twfPara.getCoeff(HF) - trial.getCoeff(HF)*network.getCoeff(HF)) < epsilon);

	// and now check the derivative (w.r. to 1 coefficient)

	// define a dummy state with 1 det and 3 times the same det as coupled (all with coeff. of 1)
	State singleDet(1);
	singleDet.coeff(0) = 1.0;
	singleDet.det(0) = HF;
	singleDet.coupledDets(0) = std::vector<detType>(0,HF);
	singleDet.coupledCoeffs(0) = std::vector<coeffType>(0,1.0);

	auto outerDerivative = std::vector<coeffType>(1,1.0);

	auto nabla = twfPara.calcNablaPars(singleDet,outerDerivative);

	// the number of parameters does not change
	assert(nabla.size() == network.pars().size());

	// and check if the derivative has the correct values
	auto baseNabla = network.calcNablaPars(singleDet,outerDerivative);

	for(size_t i = 0; i < nabla.size(); ++i){
		assert( std::abs(baseNabla(i)*trial.getCoeff(singleDet.det(0)) - nabla(i)) < epsilon);
	}

	// try something more involved
	twfPara = TrialWfPara<VecCType>(network,twfPara);
}


