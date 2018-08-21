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
#include "../Hamiltonian/ExcitationGenerators/ExcitationGenerator.hpp"
#include "../Hamiltonian/ExcitationGenerators/WeightedExcitgen.hpp"
#include "../Hamiltonian/TwoBodyHamiltonian.hpp"

namespace networkVMC{

// explicit constructor
template <typename F, typename coeffType>
Sampler<F, coeffType>::Sampler(ExcitationGenerator const &eG_,
		  detType const &HF, int numDets_):
	  excitGen(eG_.clone()),originalRef(HF),
	  numDets(numDets_),cDet(HF){};

//---------------------------------------------------------------------------//

// construct the ExcitationGenerator implicitly from the Hamiltonian
template <typename F, typename coeffType>
Sampler<F, coeffType>::Sampler(Hamiltonian const &H_,
		  detType const &HF, int numDets_):
			  excitGen(getDefaultExcitgen(H_,HF).release()),
			  originalRef(HF),numDets(numDets_),
			  cDet(HF){};

//---------------------------------------------------------------------------//

// Just take any connected determinant and let the Excitation generator decide what is connected
template <typename F, typename coeffType>
detType Sampler<F, coeffType>::getRandomConnection(detType const &startingPoint, double &p) const{
	return excitGen->generateExcitation(startingPoint, p);
}

//---------------------------------------------------------------------------//

template <typename F, typename coeffType>
double Sampler<F, coeffType>::getConnectionProb(detType const &source, detType const &target) const{
	return excitGen->getExcitationProb(source,target);
}

//---------------------------------------------------------------------------//

// Go through a list of dets and remove duplicates (generic function)
void removeDuplicate(std::vector<detType> &list){
 std::sort( list.begin(), list.end() );
// make use of std::unique for a sorted list
 auto it = std::unique( list.begin(), list.end() );
 list.erase( it, list.end() );
}

//---------------------------------------------------------------------------//

// This creates some random determinant from a given Basis
detType getRandomDeterminant(Basis const &fullBasis){
  // Essentially, we just get a random number and then return the corresponding det
  std::random_device rd;     // only used once to initialise (seed) engine
  std::mt19937 rng{rd()};    // random-number engine used (Mersenne-Twister in this case)
  std::uniform_int_distribution<int> uni(0,fullBasis.size()-1); // guaranteed unbiased
  auto randomDet = fullBasis.getDetByIndex(uni(rng));
  return randomDet; 
}

}
