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

using namespace networkVMC;

void runMetropolisTest(Sampler const &tst){
	// test the sampler tst by generating a bunch of
	// dets and checking their coeffs
	  coeffType cCoeff = coeffType();
	  coeffType lastCoeff = coeffType();
	  detType cDet = detType();
	  for(int i = 0; i<100; ++i){
		  lastCoeff = cCoeff;
		  tst.iterate(cCoeff,cDet,i);
		  std::cout<<cCoeff<<" with fraction "<<lastCoeff/cCoeff<<std::endl;
	  }
}

int main(){
  // set up the system
  int numSites{6};
  // then the Hamiltonian
  Basis basis = generateDefaultBasis(numSites);
  auto modelHam = generateDefaultHubbard(numSites);
  double U{4.}, t{-1};
  // The cost function does not matter, we only need the para to get
  // the coeffs
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


