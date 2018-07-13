/*
 * MarkovSampler.cxx
 *
 *  Created on: Jan 23, 2018
 *      Author: guther
 */

#include "MetropolisSampler.hpp"

#include <random>
#include <cmath>
#include <iostream>

#include "../Network/Parametrization.hpp"

namespace networkVMC{

template<typename T>
MetropolisSampler<T>::~MetropolisSampler() {
}

//---------------------------------------------------------------------------//

template<typename T>
void MetropolisSampler<T>::iterate(coeffType &cI, detType &dI, double &weight, 
    int i) const{
	// assign the output from the previous iteration
	cI = lastCoeff;
	dI = cDet;
  weight = 1;
	// and then get the one for the next one - this way we ensure the first iterate() call
	// returns the starting point
	// set up the rng
	std::random_device rd;     // only used once to initialise (seed) engine
	std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
	std::uniform_real_distribution<double> uni;		// Uniform distribution from 0.0 to 1.0
	// forward- and backwards generation probabilities
	double pEx, pBack;
	// First, get a random coupled determinant (from cDet)
  // need to unbias for excitgen
	// And its coefficient
  // std::norm returns |a+ib|^2
  //double prob = std::norm(tmpCoeff)/std::norm(lastCoeff);
	detType tmp{getRandomConnection(cDet,pEx)};
	// And its coefficient
	coeffType tmpCoeff{para->getCoeff(tmp)};
	// unbiasing with generation probability in principle necessary (unless prob. is symmetric)
	pBack = getConnectionProb(tmp,cDet);
	if(uni(rng) < (pBack*std::norm(tmpCoeff))/(pEx*std::norm(lastCoeff))){
	//if(uni(rng) < (std::norm(tmpCoeff/lastCoeff))){
		// With probability |cJ|^2/|cI|^2, accept the move
		cDet = tmp;
		lastCoeff = tmpCoeff;
		// and set the output
	}
}
//instantiate class
template class MetropolisSampler<VecType>;
template class MetropolisSampler<VecCType>;
}
