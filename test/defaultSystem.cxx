/*
 * defaultSystem.cxx
 *
 *  Created on: May 4, 2018
 *      Author: guther
 */

#include "../src/NNWLib.hpp"

using namespace networkVMC;

FermiHubbardHamiltonian generateDefaultHubbard(int numSites){
 // Hamiltonian setup
  int numStates = 2*numSites;
  double U{4.}, t{-1};
  return generateFermiHubbard(numStates, U, t);
}

FermionBasis generateDefaultBasis(int numSites){
  int numStates = 2*numSites;
  int spinUp = numSites/2;
  int spinDown = numSites/2;
  SpinConfig spinConfig(spinUp, spinDown, numStates);
  FermionBasis basis(spinConfig);
  return basis;
}


