/*
 * Sampler.cxx
 * Copyright (c) 2017
 * Author: Ke Liao and Kai Guther
 * All rights reserved
 */
#include "Sampler.hpp"

#include <vector>
#include <random>
#include <algorithm>
#include <Eigen/Dense>
#include <iostream>

#include "../HilbertSpace/Basis.hpp"
#include "../Hamiltonian/Hamiltonian.hpp"

namespace networkVMC{

// Just take any connected determinant and let the Hamiltonian decide what is connected
detType Sampler::getRandomConnection(detType const &startingPoint) const{
	double p{0};
	return H->getRandomCoupledState(startingPoint, p);
}


// Go through a list of dets and remove duplicates (generic function)
void removeDuplicate(std::vector<detType> &list){
 std::sort( list.begin(), list.end() );
// make use of std::unique for a sorted list
 auto it = std::unique( list.begin(), list.end() );
 list.erase( it, list.end() );
}

// This creates some random determinant from a given Basis
detType getRandomDeterminant(Basis const &fullBasis){
  // Essentially, we just get a random number and then return the corresponding det
  std::random_device rd;     // only used once to initialise (seed) engine
  std::mt19937 rng{rd()};    // random-number engine used (Mersenne-Twister in this case)
  std::uniform_int_distribution<int> uni(0,fullBasis.getSize()); // guaranteed unbiased
  auto randomDet = fullBasis.getDetByIndex(uni(rng));
  return randomDet; 
}

}
