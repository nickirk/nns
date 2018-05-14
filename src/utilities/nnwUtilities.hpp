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
#include "State.hpp"

namespace networkVMC{

// Try to match the neural network to a given input determinant

void preTrain(Parametrization<VecType> &network, State const &target, Sampler const &msampler,double trainRate=0.1);

}

#endif /* SRC_UTILITIES_NNWUTILITIES_HPP_ */
