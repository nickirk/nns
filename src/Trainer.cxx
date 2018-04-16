/*
 * Trainer.cxx
 *
 *  Created on: Jan 23, 2018
 *      Author: guther
 */

#include "Trainer.hpp"
#include <iostream>

Trainer::Trainer(Hamiltonian const &modelHam_, NeuralNetwork &NNW_, Sampler  &msampler_):modelHam(modelHam_),NNW(NNW_), msampler(msampler_) {
	int numDets{msampler.getNumDets()};
	inputState.resize(numDets);

	//optional
	//coupledDets.resize(numDets);
	//coupledCoeffs.resize(numDets);
}

//---------------------------------------------------------------------------------------------------//

Trainer::~Trainer() {
}

//---------------------------------------------------------------------------------------------------//

void Trainer::train(double learningRate, int method, int iteration){
	//int const method{2};
	// Get the number of samples
	int numDets{msampler.getNumDets()};
	inputState.clear();
	inputState.resize(numDets);
        //coupledCoeffsEpoch.clear();
        //coupledDetsEpoch.clear();
	//sampledDet.resize(numDets);
	//sampledCoeff.resize(numDets);
	//Get the first coefficient + determinant
	//inputState[0].det = msampler.getDet();
	//inputState[0].coeff = NNW.getCoeff(inputState[0].det);

	// And now, for the chosen number of samples, get the respective determinants and
	// coefficients
	for(int i=0; i < numDets; ++i){
	  msampler.iterate(inputState[i].coeff, inputState[i].det);
	  //-------------------------------------
    inputState[i].coupledDets = modelHam.getCoupledStates(inputState[i].det);
	  inputState[i].coupledCoeffs.resize(inputState[i].coupledDets.size());
	  //sampledDets[i] =  msampler.getDet(i);
	  //sampledCoeffs[i] =  NNW.getCoeff(sampledDets[i]);
          //NNW.cacheNetworkState();
	  for(size_t j=0; j < inputState[i].coupledDets.size(); ++j){
	  	inputState[i].coupledCoeffs[j]=NNW.getCoeff(inputState[i].coupledDets[j]);
	  }
	}
	// sort the list of determinants so that we can avoid 
        // computing repeated tasks.
        //coupledCoeffsEpoch.push_back(coupledCoeffs);
        //coupledDetsEpoch.push_back(coupledDets);
        //msampler.setReference(sampledDets[numDets-1]);
	std::sort(inputState.begin(), inputState.end());
	NNW.updateParameters(method,inputState,learningRate,iteration);
}
//--------------------------------------------------------------------------------------------------//
/*
std::vector<coeffType > NeuralNetwork::getCoupledCoeffs(detType const &det,
		std::vector<detType > &coupledDets) const{
	coupledDets = getCoupledStates(det);
	std::vector<coeffType > coupledCoeffs(coupledDets.size());
	for(size_t i=0; coupledDets.size(); ++i){
		feedForward(coupledDets[i]);
		coupledCoeffs[i] = outputLayer();
	}
	return coupledCoeffs;
}
*/
//---------------------------------------------------------------------------------------------------//

double Trainer::getE() const{
	// Here, we just output the value of the cost function (usually the energy) of the
	// network
	return NNW.getCostFunction()->calc(inputState);
}

//---------------------------------------------------------------------------------------------------//

std::vector<State > Trainer::getState() const{
	// This is for testing and debugging purpose, only
	return inputState;
}

