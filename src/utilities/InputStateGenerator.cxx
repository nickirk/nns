/*
 * inputStateGenerator.cxx
 *
 *  Created on: Jun 28, 2018
 *      Author: guther
 */

#include "InputStateGenerator.hpp"

#include "../Samplers/Sampler.hpp"
#include "../Network/Parametrization.hpp"
#include "../Hamiltonian/ExcitationGenerators/ConnectionGenerators/ConnectionGenerator.hpp"
#include <iostream>

#include "../Hamiltonian/TwoBodyHamiltonian.hpp"
namespace networkVMC {

template<typename T>
InputStateGenerator<T>::InputStateGenerator(Sampler &msampler_, Hamiltonian const &H_, Parametrization<T> const &para_):
	msampler(msampler_), H(H_), para(para_){
}

template<typename T>
InputStateGenerator<T>::~InputStateGenerator() {
}

template<typename T>
State InputStateGenerator<T>::generate(int numCons) const{
	// set up a state of the matching size
	int numDets{msampler.getNumDets()};
	// reset the sampler
	msampler.reset();
	State outputState(numDets);

#pragma omp parallel
  {
	// sampling is not threadsafe, so each thread creates it's own sampler
    thread_local std::unique_ptr<Sampler> samplerThread(msampler.clone());
#pragma omp for
	for(int i=0; i < numDets; ++i){
      // iterate the sampler: This also requires the iteration as an input, as
	  // some samplers pre-fetch the ensemble of determinants
	  samplerThread->iterate(outputState.coeff(i), outputState.det(i),outputState.weight(i), i);
	  if(numCons > 0)
	     // get some coupled determinants and their coefficients to use in the
	     // energy estimator
	     // Only required if the CF needs it
	     //outputState.coupledDets(i) = H.getCoupledStates(outputState.det(i));
		 outputState.coupledDets(i) = sampleConnections(H,outputState.det(i),numCons,outputState.coupledWeights(i));
	  else if(numCons < 0){
		  // get all the coupled states and assign weight to 1.
		  outputState.coupledDets(i) = H.getCoupledStates(outputState.det(i));
		  outputState.coupledWeights(i).resize(outputState.coupledDets(i).size());
		  for(size_t j=0; j < outputState.coupledDets(i).size(); ++j){
			outputState.coupledWeights(i)[j]=1.0/static_cast<double>(
					outputState.coupledDets(i).size());
		  }
	  }

	  // just get the coefficients from the NNW
	  outputState.coupledCoeffs(i).resize(outputState.coupledDets(i).size());
	  for(size_t j=0; j < outputState.coupledDets(i).size(); ++j){
		outputState.coupledCoeffs(i)[j]=para.getCoeff(outputState.coupledDets(i)[j]);
	  }
	 }
  }
  return outputState;
}

template class InputStateGenerator<VecType>;
template class InputStateGenerator<VecCType>;

} /* namespace networkVMC */