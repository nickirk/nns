/*
 * nnwUtilities.cxx
 *
 *  Created on: Jan 23, 2018
 *      Author: guther
 */

#include "State.hpp"
#include "CostFunction.hpp"
#include "NormCF.hpp"
#include "Trainer.hpp"
#include "Sampler.hpp"
#include <iostream>

void preTrain(NeuralNetwork &network, std::<State> const &target, Sampler const &msampler, double trainRate, int iteration){
// Trains the network to represent some state target
// We first backup the current cost function
	CostFunction const *backupCF = network.getCostFunction();
// Then, set the cost function to the L2-distance to target
	NormCF stateDistance(target);
	network.setCostFunction(stateDistance);
// Set up an initial list of determinants to train
// Caveat: All determinants not present in the state do not matter
// i.e. their coefficients are treated as unknown
	//std::vector<detType > list = target.getDets();

	// set up the trainer
	Trainer ev(network,msampler);
// Train the network
	int const maxTrainCount = 1000;
	for(int i = 0; i < maxTrainCount;++i){
		ev.train(trainRate);
		std::cout << "Distance/t" << ev.getE() << std::endl;
	}
	network.setCostFunction(*backupCF);
}

