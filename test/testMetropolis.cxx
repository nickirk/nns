/*
 * testMetropolis.cxx
 *
 *  Created on: Apr 30, 2018
 *      Author: guther
 */

#include "../src/NNWLib.hpp"
#include <vector>
#include <iostream>
#include <fstream>

using namespace networkVMC;

int main(){
  // set up the system
  int numSites(4);
  int numStates(2*numSites);
  int spinUp(2);
  int spinDown(2);
  SpinConfig spinConfig(spinUp, spinDown, numStates);
  // then the Hamiltonian
  Basis basis(spinConfig);
  FermionicHamiltonian modelHam(numStates);
  double U{4.}, t{-1};
  modelHam = generateFermiHubbard(numStates, U, t);
  // The cost function does not matter, we only need the NNW to get
  EnergyEsMarkov eCF(modelHam);
  int numHidden = 10;
  // set up the NNW
  RBM rbm(numStates, numHidden);
  detType HF=basis.getDetByIndex(0);
  // and the sampler
  MetropolisSampler<VecCType> sampler(modelHam, basis, HF,rbm);
  sampler.setNumDets(5000);

  //sampler.diffuse(list,spinConfig);
  //Setup the trainer
  double energy{0.0};
  double trainRate(0.0005);
  //AcceleratedGradientDescent<VecCType> sl(trainRate);
  ADAM<VecCType> sl(trainRate);
  Trainer<VecCType> ev(rbm, sampler, sl, eCF);
  coeffType cCoeff = coeffType();
  coeffType lastCoeff = coeffType();
  detType cDet = detType();
  std::ofstream myfile1;
  myfile1.open ("en1");
  for(int i = 0; i<5000; ++i){
	  lastCoeff = cCoeff;
    trainRate*=0.999;
    std::cout << "Training rate=" << trainRate << std::endl;
	  ev.train(trainRate);
	  energy = ev.getE();
    auto states=ev.getState();

    //for(size_t s=0; s<states.size(); ++s){
    //  std::cout << "C_" << s << "= " << states.coeff(s) << std::endl;
    //}
    std::cout << "iteration=" << i << std::endl;
    std::cout<<"Energy "<<energy<<std::endl;
    myfile1 << i << "  " << energy << std::endl;
  }
}


