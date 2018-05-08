/*
 * Trainer.cxx
 *
 *  Created on: Jan 23, 2018
 *      Author: guther
 */

#include "Trainer.hpp"
#include "Hamiltonian/Hamiltonian.hpp"
#include "Samplers/Sampler.hpp"
#include "Solvers/Solver.hpp"
#include "CostFunctions/CostFunction.hpp"
#include "Network/Parametrization.hpp"
#include <iostream>

namespace networkVMC{
template <typename T>
Trainer<T>::Trainer(Parametrization &NNW_, Sampler const &msampler_, Solver<T> &sl_, CostFunction const &cf_):
		modelHam(msampler_.getH()),NNW(NNW_), msampler(msampler_),sl(sl_),cf(cf_) {
	inputState.resize(msampler.getNumDets());
}

//---------------------------------------------------------------------------------------------------//
template <typename T>
Trainer<T>::~Trainer() {
}

//---------------------------------------------------------------------------------------------------//

template <typename T>
void Trainer<T>::train(double learningRate){
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
template <typename T>
void Trainer<T>::train(){
	int numDets{msampler.getNumDets()};
	inputState.resize(numDets);
	//coupledCoeffsEpoch.clear();
        //coupledDetsEpoch.clear();
	//sampledDet.resize(numDets);
	//sampledCoeff.resize(numDets);
	//Get the first coefficient + determinant
	//inputState[0].det = msampler.getDet();
	//inputState[0].coeff = NNW.getCoeff(inputState[0].det);

	// And now, for the chosen number of samples, get the respective determinants and
	// coefficients
	std::cout<<"numdets "<<numDets<<'\n';
	// TODO this should not be part of the trainer, move it to somewhere in the state
	for(int i=0; i < numDets; ++i){
	  msampler.iterate(inputState.coeff(i), inputState.det(i));
	  // get some coupled determinants and their coefficients to use in the
	  // energy estimator
      inputState.coupledDets(i) = modelHam.getCoupledStates(inputState.det(i));
	  inputState.coupledCoeffs(i).resize(inputState.coupledDets(i).size());

	  // just get the coefficients from the NNW
	  for(size_t j=0; j < inputState.coupledDets(i).size(); ++j){
	  	inputState.coupledCoeffs(i)[j]=NNW.getCoeff(inputState.coupledDets(i)[j]);
	  }
	}
	updateParameters(inputState);
}

//---------------------------------------------------------------------------------------------------//

template <typename T>
void Trainer<T>::updateParameters(State const &input){
	// first, get the derivative of the cost function with respect to the
	// wavefunction coefficients
	auto dEdC = cf.nabla(input);
	// add the inner derivative to get the full derivative
	// of the cost function with respect to the parameters
	auto dEdPars = NNW.calcNablaPars(input,dEdC);
	// feed these to the solver
	sl.update(NNW.pars(),dEdPars,input);
}

//---------------------------------------------------------------------------------------------------//

template <typename T>
double Trainer<T>::getE() const{
	// Here, we just output the value of the cost function (usually the energy) of the
	// network
	return cf.calc(inputState);
}

//---------------------------------------------------------------------------------------------------//

template <typename T>
State  Trainer<T>::getState() const{
	// This is for testing and debugging purpose, only
	return inputState;
}
//instantiate class
template class Trainer<VecType>;

}
