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
#include "Network/Parametrization.hpp"
#include "CostFunctions/CostFunction.hpp"
#include "Samplers/Sampler.hpp"
#include <Eigen/Dense>
#include <memory>

namespace networkVMC{

// forward declaration for more efficient compilation

//class CostFunction;
//class Sampler;
//class Parametrization;
class Hamiltonian;

// wrapper class for optimizing parameters
// We take a parametrization, a sampler, a solver and a cost function
// then some magic happens and the parameters are optimized with respect to
// the cost function
template <typename F=std::complex<double>, typename coeffType=std::complex<double>>
class Trainer {
public:
// supply a Parametrization, a sampler, a solver and the cost function
	using T=Eigen::Matrix<F, Eigen::Dynamic, 1>;
	Trainer(Parametrization<F, coeffType> &NNW_, Sampler<F, coeffType> &msampler,
			Solver<F, coeffType> &sl_, CostFunction<F, coeffType> &cf_, Hamiltonian const &H_);
// train() tries to optimize the parameters of the NNW with respect to its cost function
	void train();
// optionally set the learning rate for this step
	void train(double learningRate);
// read-out methods for energy and coefficients
	void getNextCoeff(coeffType &cI, detType &dI);
	coeffType getE() const;
// read out the current state of the NNW
	State<coeffType> const& getState() const;
// For testing purposes: get the cost function
	CostFunction<F, coeffType> const & getCF() const {return cf;}
// update the parameters of the Parametrization
	void updateParameters(State<coeffType> const &input);
	virtual ~Trainer();
private:
	Hamiltonian const &modelHam;
	// the parameterization of the wave function
	Parametrization<F, coeffType> &NNW;
	// the sampler generating the random states (changes its state
	// when sampling)
	Sampler<F, coeffType> &msampler;
	// the solver that optimizes the parameters (changes its parameters)
	Solver<F, coeffType> &sl;
	// The cost function with respect to which we optimize (is set up to fit the sampler)
	CostFunction<F, coeffType> &cf;

	// This is only for debugging
	State<coeffType> inputState;

	// Has a reference member, so assignment is not a thing
	Trainer& operator=(Trainer const &source);
	Trainer& operator=(Trainer &&source);

};

}

#endif /* SRC_TRAINER_HPP_ */
