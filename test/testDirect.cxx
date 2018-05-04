/*
 * testDirect.cxx
 *
 *  Created on: May 2, 2018
 *      Author: guther
 */

#include <vector>
#include <iostream>
#include "../src/NNWLib.hpp"
#include "defaultSystem.hpp"

using namespace networkVMC;

int main(){

// Hamiltonian setup
  int numSites = 8;
  auto modelHam = generateDefaultHubbard(numSites);
  auto basis = generateDefaultBasis(numSites);
// Cost function setup
  EnergyEstimator eCF(modelHam);

// Solver setup
  double trainRate(0.01);
  ADAM sl(trainRate);

// Parametrization setup
  DirectParametrization par(basis);

// Sampler setup
  detType HF=basis.getDetByIndex(0);
  ListGen sample(modelHam,basis,HF,par);

// And optimize the direct parametrization
  Trainer ev(par,sample,sl,eCF);

  for(int i=0;i<10000;++i){
	  ev.train();
	  std::cout<<"Energy "<<ev.getE()<<" in iteration "<<i<<std::endl;
  }
}



