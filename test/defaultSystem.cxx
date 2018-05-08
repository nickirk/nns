/*
 * defaultSystem.cxx
 *
 *  Created on: May 4, 2018
 *      Author: guther
 */

#include "../src/NNWLib.hpp"

using namespace networkVMC;

FermionicHamiltonian generateDefaultHubbard(int numSites){
 // Hamiltonian setup
  int numStates = 2*numSites;
  FermionicHamiltonian modelHam(numStates);
  double U{4.}, t{-1};
  modelHam = generateFermiHubbard(numStates, U, t);
  return modelHam;
}

Basis generateDefaultBasis(int numSites){
  int numStates = 2*numSites;
  int spinUp = numSites/2;
  int spinDown = numSites/2;
  SpinConfig spinConfig(spinUp, spinDown, numStates);
  Basis basis(spinConfig);
  return basis;
}


