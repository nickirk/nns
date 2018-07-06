/*
 * testSubSpaceCF.cxx
 *
 *  Created on: May 4, 2018
 *      Author: guther
 */

#include "../src/NNWLib.hpp"
#include "defaultSystem.hpp"
#include <iostream>

using namespace networkVMC;

int main(){
  int numSites = 8;
  auto modelHam = generateDefaultHubbard(numSites);
  auto basis = generateDefaultBasis(numSites);

  SubspaceCF sCF(modelHam);

  // Solver setup
  double trainRate(0.01);
  ADAM<> sl(trainRate);

  // Parametrization setup
  DirectParametrization<> par(basis);

  // Sampler setup
  detType HF=basis.getDetByIndex(0);
  ListGen<> sample(modelHam,basis,HF,par);

  // And optimize the direct parametrization
  Trainer<VecType> ev(par,sample,sl,sCF);

  for(int i=0;i<1000;++i){
	  ev.train();
	  std::cout<<"Energy "<<ev.getE()<<" in iteration "<<i<<std::endl;
  }
}


