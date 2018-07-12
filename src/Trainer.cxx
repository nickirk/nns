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
#include "Network/Parametrization.hpp"
#include <iostream>
#include <memory>
#include "CostFunctions/CostFunction.hpp"
#include "utilities/InputStateGenerator.hpp"

namespace networkVMC{

// TODO: Add more constructors, with default arguments
template <typename T>
Trainer<T>::Trainer(Parametrization<T> &NNW_, Sampler const &msampler_,
		Solver<T> &sl_, CostFunction &cf_, Hamiltonian const& H_):
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
	// the sampler dictates how many determinants we use
	int numDets{msampler.getNumDets()};
	inputState.resize(numDets);
  int spaceSize =0;

	// And now, for the chosen number of samples, get the respective determinants and
	// coefficients
/*
<<<<<<< HEAD
#pragma omp parallel
	{
	// sampling is not threadsafe, so each thread creates it's own sampler
    std::unique_ptr<Sampler> samplerThread(msampler.clone());
#pragma omp for reduction(+:spaceSize)
	for(int i=0; i < numDets; ++i){
    // iterate the sampler: This also requires the iteration as an input, as
	  // some samplers pre-fetch the ensemble of determinants
	  samplerThread->iterate(
        inputState.coeff(i), inputState.det(i), inputState.weight(i),i
        );
	  // get some coupled determinants and their coefficients to use in the
	  // energy estimator
    inputState.coupledDets(i) = modelHam.getCoupledStates(inputState.det(i));
	  inputState.coupledCoeffs(i).resize(inputState.coupledDets(i).size());
	  inputState.coupledWeights(i).resize(inputState.coupledDets(i).size());
    spaceSize+=1;
	  // just get the coefficients from the NNW
    // also set the weight of the coupled dets
	  for(size_t j=0; j < inputState.coupledDets(i).size(); ++j){
	  	inputState.coupledCoeffs(i)[j]=NNW.getCoeff(inputState.coupledDets(i)[j]);
        // at this stage every weight is 1.
      inputState.coupledWeights(i)[j]=1;
      spaceSize+=1;
	  }
	}
	}
  inputState.spaceSize = spaceSize;
  std::cout << "spaceSize=" << spaceSize << std::endl;
  //inputState.reduce();
  // count the repeatations, inside State? 
  // count the repeated determinants, update the weights or use hash table
  // in the sampler to count the number. TODO
=======
*/
	std::cout<<"numdets "<<numDets<<'\n';

	// create the input State
	InputStateGenerator<T> isg(msampler,modelHam, NNW);
	// the number of connections required is passed via the cost function
	std::cout << "Trainer.cxx: numCons=" << cf.connectionsRequired() << std::endl;
	inputState = isg.generate(cf.connectionsRequired());

	// use the data obtained in generation of the input state to
	// set new biases for the sampling
	msampler.updateBiases();

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
	T dEdPars;
	switch(msampler.type()){
	case Markov:
		// for markov-type samplers, use EnergyEsMarkov
		dEdPars = NNW.calcNablaParsConnected(input,dEdC);
		break;
	case PreFetched:
		dEdPars = NNW.calcNablaPars(input,dEdC);
		break;
	default:
		 dEdPars = NNW.calcNablaParsConnected(input,dEdC);
	}
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
State const &Trainer<T>::getState() const{
	// This is for testing and debugging purpose, only
	return inputState;
}
//instantiate class
template class Trainer<VecType>;

//instantiate class
template class Trainer<VecCType>;
}
