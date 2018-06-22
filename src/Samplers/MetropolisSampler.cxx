/*
 * MarkovSampler.cxx
 *
 *  Created on: Jan 23, 2018
 *      Author: guther
 */

#include "MetropolisSampler.hpp"

#include <random>
#include <cmath>

#include "../Network/Parametrization.hpp"

namespace networkVMC{

template<typename T>
MetropolisSampler<T>::~MetropolisSampler() {
}

template<typename T>
void MetropolisSampler<T>::iterate(coeffType &cI, detType &dI, double &weight, 
    int i) const{
	// set up the rng
	std::random_device rd;     // only used once to initialise (seed) engine
	std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
	std::uniform_real_distribution<double> uni;		// Uniform distribution from 0.0 to 1.0
	// First, get a random coupled determinant (from cDet)
  // need to unbias for excitgen
	detType tmp{getRandomConnection(cDet)};
	// And its coefficient
	coeffType tmpCoeff{para->getCoeff(tmp)};
  // std::norm returns |a+ib|^2
  double prob = std::norm(tmpCoeff)/std::norm(lastCoeff);
	if(uni(rng) < prob){
		// With probability cJ/cI, accept the move
		cDet = tmp;
		lastCoeff = tmpCoeff;
		// and set the output
	}
	cI = lastCoeff;
	dI = cDet;
  weight = 1;
}
//instantiate class
template class MetropolisSampler<VecType>;
template class MetropolisSampler<VecCType>;
}
