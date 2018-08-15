/*
 * MarkovSampler.cxx
 *
 *  Created on: Jan 23, 2018
 *      Author: guther and Liao
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
  // double prob = std::norm(tmpCoeff)/std::norm(lastCoeff);
  detType tmp;
  double prob(0.);
	coeffType tmpCoeff;
  double probRand =std::max(100*std::pow(0.9,i), 0.);
  if (rng() < probRand) {
    tmp = getRandomDeterminant(*fullBasis);
    //tmp=getRandomConnection(cDet,pEx);
    // set pEx and pBack both to 1, because the back and forth prob are the same
    // and only the ratio matters.
	  pEx = 1.; 
    pBack = 1.;
  }
  else {
    tmp=getRandomConnection(cDet,pEx);
	  pBack = getConnectionProb(tmp,cDet);
  }
  tmpCoeff = para->getCoeff(tmp);
  //coeffType tmpCoeffDeref;
  //coeffType lastCoeffDeref;
  //if (verbatimCast(tmp)==15) tmpCoeffDeref=1.0;
  //else tmpCoeffDeref=std::pow(0.5,2);
  //if (verbatimCast(cDet)==15) lastCoeffDeref=1.0;
  //else lastCoeffDeref=std::(0.5,2);
	// And its coefficient
	// unbiasing with generation probability in principle necessary (unless prob. is symmetric)

  //prob = pBack/pEx*std::norm(tmpCoeffDeref/lastCoeffDeref);
  prob = pBack/pEx*std::norm(tmpCoeff/lastCoeff);
	if(rng() < prob){
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
