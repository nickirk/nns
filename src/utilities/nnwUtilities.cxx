/*
 * nnwUtilities.cxx
 *
 *  Created on: Jan 23, 2018
 *      Author: guther
 */

#include "../Trainer.hpp"
#include <iostream>

#include "../CostFunctions/CostFunction.hpp"
#include "../CostFunctions/NormCF.hpp"
#include "../Samplers/Sampler.hpp"
#include "../Solvers/ADAM.hpp"
#include "State.hpp"

namespace networkVMC{

void preTrain(Parametrization &network, State const &target, Sampler const &msampler, double trainRate){
// Trains the network to represent some state target
// Then, set the cost function to the L2-distance to target
	NormCF stateDistance(target);
// Set up an initial list of determinants to train
// Caveat: All determinants not present in the state do not matter
// i.e. their coefficients are treated as unknown
	//std::vector<detType > list = target.getDets();
	ADAM sl(trainRate);
	// set up the trainer
	Trainer ev(network,msampler,sl,stateDistance);
    // Train the network
	int const maxTrainCount = 1000;
	for(int i = 0; i < maxTrainCount;++i){
		ev.train();
		std::cout << "Distance/t" << ev.getE() << std::endl;
	}
}

}

