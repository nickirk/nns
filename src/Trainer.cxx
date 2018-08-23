/*
 * Trainer.cxx
 *
 *  Created on: Jan 23, 2018
 *      Author: guther
 */

#include "Trainer.hpp"
#include "Samplers/Sampler.hpp"
#include "Solvers/Solver.hpp"
#include "Network/Parametrization.hpp"
#include <iostream>
#include <memory>
#include "CostFunctions/CostFunction.hpp"
#include "Hamiltonian/TwoBodyHamiltonian.hpp"
#include "utilities/InputStateGenerator.hpp"

namespace networkVMC{

// TODO: Add more constructors, with default arguments
template <typename F, typename coeffType>
Trainer<F, coeffType>::Trainer(Parametrization<F, coeffType> &NNW_, Sampler<coeffType> &msampler_,
		Solver<F, coeffType> &sl_, CostFunction<F, coeffType> &cf_, Hamiltonian const& H_):
		modelHam(H_),NNW(NNW_), msampler(msampler_),sl(sl_),
		// here, we make sure that the cost function and the
		// sampler are compatible
		cf(cf_) {
	// make sure the cost function and the sampler are compatible
	cf.setUpCF(msampler.type());
	// prepare the input state
	inputState.resize(msampler.getNumDets());
}

//---------------------------------------------------------------------------------------------------//
template <typename F, typename coeffType>
Trainer<F, coeffType>::~Trainer() {
}

//---------------------------------------------------------------------------------------------------//

template <typename F, typename coeffType>
void Trainer<F, coeffType>::train(double learningRate){
	// only change the learningRate for this iteration
	auto tmp = sl.getLearningRate();
	sl.setLearningRate(learningRate);
	// do the normal train()
	train();
	// then reset the learningRate
	sl.setLearningRate(tmp);
}

//---------------------------------------------------------------------------------------------------//


// prepare an input
template <typename F, typename coeffType>
void Trainer<F, coeffType>::train(){
	// the sampler dictates how many determinants we use
	int numDets = msampler.getNumDets();
	inputState.resize(numDets);

	std::cout<< "Trainer.cxx: numdets= " << numDets << std::endl;

	// create the input State
	InputStateGenerator<F, coeffType> isg(msampler,modelHam, NNW);
	// the number of connections required is passed via the cost function
	inputState = isg.generate(cf.connectionsRequired());

	// use the data obtained in generation of the input state to
	// set new biases for the sampling
	msampler.updateBiases();

	updateParameters(inputState);
}

//---------------------------------------------------------------------------------------------------//

template <typename F, typename coeffType>
void Trainer<F, coeffType>::updateParameters(State<coeffType> const &input){
	// first, get the derivative of the cost function with respect to the
	// wavefunction coefficients
	auto dEdC = cf.nabla(input);
	// add the inner derivative to get the full derivative
	// of the cost function with respect to the parameters
	T dEdPars;
	switch(msampler.type()){
	case Markov:
		// for markov-type samplers, use EnergyEsMarkov
		dEdPars = NNW.calcNablaParsMarkovConnected(input,dEdC,getE());
		break;
	case PreFetched:
		dEdPars = NNW.calcNablaParsConnected(input,dEdC);
		break;
	default:
		throw SamplerTypeDoesNotExist(msampler.type());
	}
	// feed these to the solver
	sl.update(NNW.pars(),dEdPars,input,msampler.type());
}

//---------------------------------------------------------------------------------------------------//

template <typename F, typename coeffType>
coeffType Trainer<F, coeffType>::getE() const{
	// Here, we just output the value of the cost function (usually the energy) of the
	// network
	return cf.calc(inputState);
}

//---------------------------------------------------------------------------------------------------//

template <typename F, typename coeffType>
State<coeffType> const &Trainer<F, coeffType>::getState() const{
	// This is for testing and debugging purpose, only
	return inputState;
}
//instantiate class
template class Trainer<double, double>;

//instantiate class
template class Trainer<std::complex<double>, std::complex<double>>;
}
