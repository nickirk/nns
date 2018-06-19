/*
 * testMetropolis.cxx
 *
 *  Created on: Apr 30, 2018
 *      Author: guther
 */

#include "../src/NNWLib.hpp"
#include <vector>
#include <iostream>

using namespace networkVMC;

int main(){
  // set up the system
  int numSites(6);
  int numStates(2*numSites);
  int spinUp(3);
  int spinDown(3);
  SpinConfig spinConfig(spinUp, spinDown, numStates);
  // then the Hamiltonian
  Basis basis(spinConfig);
  FermionicHamiltonian modelHam(numStates);
  double U{4.}, t{-1};
  modelHam = generateFermiHubbard(numStates, U, t);
  // The cost function does not matter, we only need the NNW to get
  EnergyEstimator eCF(modelHam);
  NeuralNetwork<> NNW(eCF);
  // set up the NNW
  NNW.constrInputLayer(numStates);
  NNW.constrDenseLayer(NNW.getLayer(0)->getActs(), "Tanh", 10*numStates);
  NNW.constrDenseLayer(NNW.getLayer(1)->getActs(), "Linear", 2);
  NNW.initialiseNetwork();

  detType HF=basis.getDetByIndex(0);
  // and the sampler
  MetropolisSampler<> sampler(modelHam, basis, HF,NNW);

  coeffType cCoeff = coeffType();
  coeffType lastCoeff = coeffType();
  detType cDet = detType();
  for(int i = 0; i<100; ++i){
	  lastCoeff = cCoeff;
	  sampler.iterate(cCoeff,cDet,i);
	  std::cout<<cCoeff<<" with fraction "<<lastCoeff/cCoeff<<std::endl;
  }
}


