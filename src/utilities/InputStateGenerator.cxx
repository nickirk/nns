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
#include "../Hamiltonian/TwoBodyHamiltonian.hpp"
#include "../HilbertSpace/Determinant.hpp"

#include <memory>
#include <iostream>
namespace networkVMC {

InputStateGenerator::InputStateGenerator(Sampler &msampler_, Hamiltonian const &H_,
		Parametrization const &para_):
	msampler(msampler_), H(H_), para(para_){
}

InputStateGenerator::~InputStateGenerator() {
}

State InputStateGenerator::generate(int numCons) const{
  // set up a state of the matching size
  int numDets{msampler.getNumDets()};
  State outputState(numDets);
  int accept(0);
  #pragma omp parallel
  {
    // sampling is not threadsafe, so each thread creates it's own sampler
    std::unique_ptr<Sampler> samplerThread(msampler.clone());

    #pragma omp for
    for(int i=0; i < numDets; ++i){
      // iterate the sampler: This also requires the iteration as an input, as
      // some samplers pre-fetch the ensemble of determinants
      samplerThread->iterate(outputState.coeff(i), outputState.det(i),outputState.weight(i), i);
      //accumulate the totalWeights for later use
      outputState.addTotalWeights(outputState.weight(i));
      if (i==0) accept=0;
      else if (outputState.det(i-1) != outputState.det(i)) accept++;
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
  std::cout << "InputStateGenerator.cxx: Accpt rate=" << double(accept)/numDets << std::endl;
  return outputState;
}

} /* namespace networkVMC */
