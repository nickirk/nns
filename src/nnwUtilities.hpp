/*
 * nnwUtilities.hpp
 *
 *  Created on: Jan 23, 2018
 *      Author: guther
 */

#ifndef SRC_NNWUTILITIES_HPP_
#define SRC_NNWUTILITIES_HPP_

#include "Nnw.hpp"
#include "Sampler.hpp"
#include "State.hpp"

// Try to match the neural network to a given input determinant

void preTrain(NeuralNetwork &network, std::<State> const &target, Sampler const &msampler, int iteration);

#endif /* SRC_NNWUTILITIES_HPP_ */
