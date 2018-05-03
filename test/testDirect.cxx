/*
 * testDirect.cxx
 *
 *  Created on: May 2, 2018
 *      Author: guther
 */

#include <vector>
#include <iostream>
#include "../src/NNWLib.hpp"

using namespace networkVMC;

int main(){

// Hamiltonian setup
  int numSites(6);
  int numStates(2*numSites);
  int spinUp(3);
  int spinDown(3);
  SpinConfig spinConfig(spinUp, spinDown, numStates);
  Basis basis(spinConfig);
  FermionicHamiltonian modelHam(numStates);
  double U{4.}, t{-1};
  modelHam = generateFermiHubbard(numStates, U, t);

// Cost function setup
  EnergyCF eCF(modelHam);

// Solver setup
  double trainRate(0.1);
  ADAM sl(trainRate);

// Parametrization setup
  DirectParametrization par(basis);

// Sampler setup
  detType HF=basis.getDetByIndex(0);
  ListGen sample(modelHam,basis,HF,par,1000);

// And optimize the direct parametrization
  Trainer ev(par,sample,sl,eCF);

  for(int i=0;i<1000;++i){
	  ev.train();
	  std::cout<<"Energy "<<ev.getE()<<std::endl;
  }
}



