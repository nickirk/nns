/*
 * nnwUtilities.hpp
 *
 *  Created on: Jan 23, 2018
 *      Author: guther
 */

#ifndef SRC_UTILITIES_NNWUTILITIES_HPP_
#define SRC_UTILITIES_NNWUTILITIES_HPP_

#include "../Network/Nnw.hpp"
#include "../Samplers/Sampler.hpp"
#include "../CostFunctions/NormCF.hpp"
#include "../Solvers/ADAM.hpp"
#include "../Trainer.hpp"
#include "../Hamiltonian/FermiHubbardHamiltonian.hpp"
#include "StateForward.hpp"
#include <iostream>

namespace networkVMC{

// Try to match the neural network to a given input determinant

void preTrain(Parametrization &network, State const &target, Sampler const &msampler, double trainRate){
// Trains the network to represent some state target
// Then, set the cost function to the L2-distance to target
	NormCF stateDistance(target);
	// dummy hamiltonian, we need to get rid of this TODO
	FermiHubbardHamiltonian dummy(1.0,1.0,1);
// Set up an initial list of determinants to train
// Caveat: All determinants not present in the state do not matter
// i.e. their coefficients are treated as unknown
	//std::vector<detType > list = target.getDets();
	ADAM sl(trainRate);
	// set up the trainer
	Trainer ev(network,msampler,sl,stateDistance,dummy);
    // Train the network
	int const maxTrainCount = 1000;
	for(int i = 0; i < maxTrainCount;++i){
		ev.train();
		std::cout << "Distance/t" << ev.getE() << std::endl;
	}
}

}

#endif /* SRC_UTILITIES_NNWUTILITIES_HPP_ */
