/*
 * testPreTrain.cxx
 *
 *  Created on: Oct 26, 2017
 *      Author: guther
 */

#include "../src/CostFunctions/EnergyCF.hpp"
#include "../src/utilities/nnwUtilities.hpp"
#include "../src/utilities/State.hpp"
#include "../src/HilbertSpace/Determinant.hpp"
#include "../src/utilities/TypeDefine.hpp"
#include "../src/Network/Nnw.hpp"
#include "../src/Hamiltonian/FermiHubbardHamiltonian.hpp"
#include "../src/HilbertSpace/BasisGenerator.hpp"
#include "../src/Samplers/MetropolisSampler.hpp"
#include "../src/utilities/SpinConfig.hpp"

using namespace networkVMC;

int main(){
	int numSites{3}, numStates{2*numSites};;
	double U{4}, t{-1};
	auto modelHam = FermiHubbardHamiltonian(U,t,numSites);
	EnergyCF eCF(modelHam);
// Two types of coefficients: 0 and 1 in the target state
	coeffType coeffZero;
	coeffType coeffOne = coeffZero;
	coeffOne = 1.0;

// Now generate the basis states for the target
	SpinConfig sC(numSites,numSites,numStates);
	Basis basisGen(sC);
	detType D(numStates,false);
	std::vector<detType > targetDets(basisGen.size(), D);
	std::vector<coeffType > targetCoeffs(basisGen.size(), coeffZero);
	targetCoeffs[0] = coeffOne;

	State targetState(targetDets, targetCoeffs);

	NeuralNetwork network(eCF);
	MetropolisSampler msampler(modelHam,basisGen,D,network);
	preTrain(network,targetState,msampler);
}


