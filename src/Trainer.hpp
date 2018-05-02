/*
 * Trainer.hpp
 *
 *  Created on: Jan 23, 2018
 *      Author: guther
 */

#ifndef SRC_TRAINER_HPP_
#define SRC_TRAINER_HPP_

#include "Network/Parametrization.hpp"
#include "HilbertSpace/Determinant.hpp"
#include "Samplers/Sampler.hpp"
#include "utilities/State.hpp"
#include "utilities/TypeDefine.hpp"
#include "CostFunctions/CostFunction.hpp"
#include "Solvers/Solver.hpp"

namespace networkVMC{

// wrapper class for optimizing parameters
// We take a parametrization, a sampler, a solver and a cost function
// then some magic happens and the parameters are optimized with respect to
// the cost function
class Trainer {
public:
// supply a sampler, a Hamiltonian and the Parametrization
	Trainer(Parametrization &NNW_, Sampler const &msampler, Solver &sl, CostFunction const &cf);
// train() tries to optimize the parameters of the NNW with respect to its cost function
	void train();
// optionally set the learning rate for this step
	void train(double learningRate);
// read-out methods for energy and coefficients
	void getNextCoeff(coeffType &cI, detType &dI);
	double getE() const;
// read out the current state of the NNW
	State getState() const;
// update the parameters of the Parametrization
	void updateParameters(State const &input);
	virtual ~Trainer();
private:
	Hamiltonian const &modelHam;
	// the parameterization of the wave function
	Parametrization &NNW;
	Sampler const &msampler;
	// the solver that optimizes the parameters (changes its parameters)
	Solver &sl;
	CostFunction const &cf;

	// This is only for debugging
	State inputState;

};

}

#endif /* SRC_TRAINER_HPP_ */
