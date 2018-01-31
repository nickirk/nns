/*
 * Trainer.cxx
 *
 *  Created on: Jan 23, 2018
 *      Author: guther
 */

#include "Trainer.hpp"
#include <iostream>

Trainer::Trainer(NeuralNetwork &NNW_, Sampler  &msampler_):NNW(NNW_), msampler(msampler_) {
	int numDets{msampler.getNumDets()};
	sampledDets.resize(numDets);
	sampledCoeffs.resize(numDets);

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
        coupledCoeffsEpoch.clear();
        coupledDetsEpoch.clear();
	sampledDets.resize(numDets);
	sampledCoeffs.resize(numDets);
	//Get the first coefficient + determinant
	sampledDets[0] = msampler.getDet();
	sampledCoeffs[0] = NNW.getCoeff(sampledDets[0]);

        std::vector<detType > coupledDets = getCoupledStates(sampledDets[0]); 
	std::vector<coeffType > coupledCoeffs(coupledDets.size());
	for(size_t i=0; i < coupledDets.size(); ++i){
	  coupledCoeffs[i]=NNW.getCoeff(coupledDets[i]);
	}
        coupledCoeffsEpoch.push_back(coupledCoeffs);
        coupledDetsEpoch.push_back(coupledDets);
	// And now, for the chosen number of samples, get the respective determinants and
	// coefficients
	for(int i=1; i < numDets; ++i){
	  msampler.iterate(sampledCoeffs[i],sampledDets[i]);
          coupledDets = getCoupledStates(sampledDets[i]); 
	  coupledCoeffs.resize(coupledDets.size());
	  //sampledDets[i] =  msampler.getDet(i);
	  //sampledCoeffs[i] =  NNW.getCoeff(sampledDets[i]);
          //NNW.cacheNetworkState();
	  for(size_t j=0; j < coupledDets.size(); ++j){
	  	coupledCoeffs[j]=NNW.getCoeff(coupledDets[j]);
	  }
          coupledCoeffsEpoch.push_back(coupledCoeffs);
          coupledDetsEpoch.push_back(coupledDets);
	}
        //msampler.setReference(sampledDets[numDets-1]);
	State inputState(sampledDets,sampledCoeffs,coupledDetsEpoch,coupledCoeffsEpoch);
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
	State inputState(sampledDets,sampledCoeffs,coupledDetsEpoch,coupledCoeffsEpoch);
	return NNW.getCostFunction()->calc(inputState);
}

//---------------------------------------------------------------------------------------------------//

State Trainer::getState() const{
	// This is for testing and debugging purpose, only
	return State(sampledDets,sampledCoeffs,coupledDetsEpoch,coupledCoeffsEpoch);
}

