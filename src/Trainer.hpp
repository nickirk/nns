/*
 * Trainer.hpp
 *
 *  Created on: Jan 23, 2018
 *      Author: guther
 */

#ifndef SRC_TRAINER_HPP_
#define SRC_TRAINER_HPP_

#include "utilities/State.hpp"
#include "utilities/TypeDefine.hpp"
#include "Solvers/Solver.hpp"

namespace networkVMC{

// forward declaration for more efficient compilation

class CostFunction;
class Sampler;
class Parametrization;
class Hamiltonian;

// wrapper class for optimizing parameters
// We take a parametrization, a sampler, a solver and a cost function
// then some magic happens and the parameters are optimized with respect to
// the cost function
template <typename T=VecType>
class Trainer {
public:
// supply a Parametrization, a sampler, a solver and the cost function
	Trainer(Parametrization &NNW_, Sampler const &msampler, Solver<T> &sl, CostFunction const &cf);
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
	Solver<T> &sl;
	CostFunction const &cf;

	// This is only for debugging
	State inputState;

	// Has a reference member, so assignment is not a thing
	Trainer& operator=(Trainer const &source);
	Trainer& operator=(Trainer &&source);

};

}

#endif /* SRC_TRAINER_HPP_ */
