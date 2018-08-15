/*
 * MarkovSampler.cxx
 *
 *  Created on: Jan 23, 2018
 *      Author: guther
 */

#include "MetropolisSampler.hpp"

#include <cmath>
#include <iostream>
#include "../utilities/RNGWrapper.hpp"
#include "../Network/Parametrization.hpp"

namespace networkVMC{

template<typename T>
MetropolisSampler<T>::~MetropolisSampler() {
}

//---------------------------------------------------------------------------//

template<typename T>
void MetropolisSampler<T>::iterate(coeffType &cI, detType &dI, double &weight, 
    int i){
  weight = 1;
	// and then get the one for the next one - this way we ensure the first iterate() call
	// returns the starting point
	// set up the rng
  	RNGWrapper rng;
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
  //std::cout << "MetropolisSampler.cxx: uni(rng)=" << p << std::endl;
  //std::cout << "MetropolisSampler.cxx: pJump=" << std::norm(tmpCoeff/lastCoeff)*pBack/pEx << std::endl;

	if(rng() < (pBack/pEx*std::norm(tmpCoeff/lastCoeff))){
	//if(p < (std::norm(tmpCoeff/lastCoeff))){
		// With probability |cJ/cI|^2, accept the move
		cDet = tmp;
		lastCoeff = tmpCoeff;
		// and set the output
	}
	// assign the output from the previous iteration
	cI = lastCoeff;
	dI = cDet;
}
//instantiate class
template class MetropolisSampler<VecType>;
template class MetropolisSampler<VecCType>;
}
