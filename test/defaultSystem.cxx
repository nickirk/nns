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
// mainly there to hide the construction for historic reasons
  double U{4.}, t{-1};
  return FermiHubbardHamiltonian(U, t, numSites);
}

SpinConfig generateDefaultSpinConfig(int numSites){
  int numStates = 2*numSites;
  int spinUp = numSites/2;
  int spinDown = numSites/2;
  return SpinConfig(spinUp, spinDown, numStates);
}

Basis generateDefaultBasis(int numSites){
	return Basis(generateDefaultSpinConfig(numSites));
}


