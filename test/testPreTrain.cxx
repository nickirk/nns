/*
 * testPreTrain.cxx
 *
 *  Created on: Oct 26, 2017
 *      Author: guther
 */

#include "../src/Nnw.hpp"
#include "../src/FermionicHamiltonian.hpp"
#include "../src/EnergyCF.hpp"
#include "../src/State.hpp"
#include "../src/Determinant.hpp"
#include "../src/CoeffType.hpp"
#include "../src/Basis.hpp"

int main(){
	int numSites{3}, numStates{2*numSites};;
	double U{4}, t{-1};
	FermionicHamiltonian modelHam = generateFermiHubbard(numStates,U,t);
	EnergyCF eCF(modelHam);
// Two types of coefficients: 0 and 1 in the target state
	coeffType coeffZero = Eigen::VectorXd::Zero(2);
	coeffType coeffOne = coeffZero;
	coeffOne[0] = 1.0;

// Now generate the basis states for the target
	Basis basisGen(numStates,numSites);
	detType D(numStates,false);
	std::vector<detType > targetDets(basisGen.getSize(), D);
	std::vector<coeffType > targetCoeffs(basisGen.getSize(), coeffZero);
	targetCoeffs[0] = coeffOne;

	State targetState(targetDets, targetCoeffs);

	std::vector<int> sizeNNW = {numStates, 2, 2};
	NeuralNetwork network(sizeNNW, eCF);
	preTrain(network, targetState, 0.0001);
}


