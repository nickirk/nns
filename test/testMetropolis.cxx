/*
 * testMetropolis.cxx
 *
 *  Created on: Apr 30, 2018
 *      Author: guther
 */

#include "../src/NNWLib.hpp"
#include "defaultSystem.hpp"
#include <vector>
#include <iostream>
#include <fstream>

//TODO : Write a statistical test that checks if the correct probability
//       density is generated

using namespace networkVMC;

void runMetropolisTest(Sampler &tst){
	// test the sampler tst by generating a bunch of
	// dets and checking their coeffs
	  coeffType cCoeff = coeffType();
	  coeffType lastCoeff = coeffType();
	  detType cDet = detType();
	  double cWeight;
	  for(int i = 0; i<100; ++i){
		  lastCoeff = cCoeff;
		  tst.iterate(cCoeff,cDet, cWeight, i);
		  std::cout<<cCoeff<<" with fraction "<<lastCoeff/cCoeff<<std::endl;
	  }
	  tst.updateBiases();
}

int main(){
  // set up the system
  int numSites{6};
  // then the Hamiltonian
  Basis basis = generateDefaultBasis(numSites);
  auto modelHam = generateDefaultHubbard(numSites);
  // The cost function does not matter, we only need the para to get the coeffs
  DirectParametrization<> para(basis);

  detType HF=basis.getDetByIndex(0);
  // and the sampler
  // run 4 tests: one for the weighted Excitgen
  WeightedExcitgen wEG(modelHam,HF);
  MetropolisSampler<> weightedSampler(wEG, basis, HF, para);

  // one for the RSHubbard
  RSHubbardExcitgen hubbardEG{};
  MetropolisSampler<> hubSampler(hubbardEG, basis, HF, para);

  runMetropolisTest(hubSampler);
  // one for RSHubbard via default initialization
  MetropolisSampler<> sampler(modelHam, basis, HF,para);
  runMetropolisTest(sampler);

  // and one for the Uniform
  UniformExcitgen uniEG(HF);
  MetropolisSampler<> ugSampler(uniEG, basis, HF, para);


  runMetropolisTest(ugSampler);
}


