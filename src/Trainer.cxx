/*
 * Trainer.cxx
 *
 *  Created on: Jan 23, 2018
 *      Author: guther
 */

#include "Trainer.hpp"

Trainer::Trainer(NeuralNetwork &NNW_, Sampler const &msampler_):NNW(NNW_), msampler(msampler_) {
	int numDets{msampler.getNumDets()};
	sampledDets.resize(numDets);
	sampledCoeffs.resize(numDets);

	//optional
	coupledDets.resize(numDets);
	coupledCoeffs.resize(numDets);
}

//---------------------------------------------------------------------------------------------------//

Trainer::~Trainer() {
}

//---------------------------------------------------------------------------------------------------//

void Trainer::train(double learningRate){
	int const method{2};
	// Get the number of samples
	int numDets{msampler.getNumDets()};

	//Get the first coefficient + determinant
	sampledDets[0] = msampler.getDet();
	sampledCoeffs[0] = NNW.getCoeff(sampledDets[0]);
	// This is a bit sloppy: We want to store the current state of the network for later
	// evaluation
	NNW.cacheNetworkState();
	// And now, for the chosen number of samples, get the respective determinants and
	// coefficients
	for(int i=1; i < numDets; ++i){
		msampler.iterate(sampledCoeffs[i],sampledDets[i]);
	}
	State inputState(sampledDets,sampledCoeffs,coupledDets,coupledCoeffs);
	NNW.updateParameters(method,inputState,learningRate);
}

//---------------------------------------------------------------------------------------------------//

double Trainer::getE() const{
	// Here, we just output the value of the cost function (usually the energy) of the
	// network
	State inputState(sampledDets,sampledCoeffs,coupledDets,coupledCoeffs);
	return NNW.getCostFunction()->calc(inputState);
}

//---------------------------------------------------------------------------------------------------//

State Trainer::getState() const{
	// This is for testing and debugging purpose, only
	return State(sampledDets,sampledCoeffs,coupledDets,coupledCoeffs);
}

