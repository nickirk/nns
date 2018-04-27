/*
 * Trainer.hpp
 *
 *  Created on: Jan 23, 2018
 *      Author: guther
 */

#ifndef SRC_TRAINER_HPP_
#define SRC_TRAINER_HPP_

#include "Network/Nnw.hpp"
#include "HilbertSpace/Determinant.hpp"
#include "Samplers/Sampler.hpp"
#include "utilities/State.hpp"
#include "utilities/TypeDefine.hpp"

namespace networkVMC{

// wrapper class for training the NNW
class Trainer {
public:
// supply a sampler, a Hamiltonian and the NNW
	Trainer(NeuralNetwork &NNW_, Sampler const &msampler);
// train() tries to optimize the parameters of the NNW with respect to its cost function
	void train(double learningRate, int method, int iteration);
// read-out methods for energy and coefficients
	void getNextCoeff(coeffType &cI, detType &dI);
	double getE() const;
// read out the current state of the NNW
	State getState() const;
	virtual ~Trainer();
private:
	Hamiltonian const &modelHam;
	NeuralNetwork &NNW;
	Sampler const &msampler;
	State inputState;

};

}

#endif /* SRC_TRAINER_HPP_ */
