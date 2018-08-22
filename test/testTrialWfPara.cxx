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
	int numStates = numSites*2;
	int nHidden = 8;
	int spinUp = 2;
	int spinDown = 2;
	AbInitioHamiltonian modelHam(0);
	std::string file_name = "FCIDUMP";
	modelHam = readAbInitioHamiltonian(file_name,1);
	numStates = modelHam.getNumOrbs();
	SpinConfig spinConfig(spinUp, spinDown,numStates);// numStates);
	Basis basis(spinConfig,modelHam);
	// do a simple test of the TrialWfPara with
	// an RBM as a base para
	RBM<> network(numStates,nHidden);
	auto HF = basis.getDetByIndex(0);

	// and another para as a trial WF
    auto trial = network;
    // but read the parameters from file
    trial.readParsFromFile("Parameters.dat");
	TrialWfPara<std::complex<double>> twfPara(network,trial);

	// first compare the coefficients of the base*trial with that of the TrialWfPara
	assert(std::abs(twfPara.getCoeff(HF) - trial.getCoeff(HF)*network.getCoeff(HF)) < epsilon);

	// and now check the derivative (w.r. to 1 coefficient)

	// define a dummy state with 1 det and 3 times the same det as coupled (all with coeff. of 1)
	State<> singleDet(1);
	singleDet.coeff(0) = 1.0;
	singleDet.det(0) = HF;
	singleDet.coupledDets(0) = std::vector<detType>(0,HF);
	singleDet.coupledCoeffs(0) = std::vector<cType>(0,1.0);

	Eigen::VectorXcd outerDerivative(1);
	outerDerivative << cType(1., 0.);

	auto nabla = twfPara.calcNablaParsConnected(singleDet,outerDerivative);

	// the number of parameters does not change
	assert(nabla.size() == network.pars().size());

	// and check if the derivative has the correct values
	auto baseNabla = network.calcNablaParsConnected(singleDet,outerDerivative);

	for(size_t i = 0; i < nabla.size(); ++i){
		assert( std::abs(baseNabla(i)/trial.getCoeff(singleDet.det(0)) - nabla(i)) < epsilon);
	}

	// try something more involved
	twfPara = TrialWfPara<std::complex<double>>(network,twfPara);
}


