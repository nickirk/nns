/*
 * MarkovSampler.cxx
 *
 *  Created on: Jan 23, 2018
 *      Author: guther and Liao
 */

#include "MetropolisSampler.hpp"

#include <cmath>
#include "../utilities/RNGWrapper.hpp"
#include "../Network/Parametrization.hpp"
#include "../HilbertSpace/Basis.hpp"
#include "../HilbertSpace/Determinant.hpp"

namespace networkVMC{

MetropolisSampler::MetropolisSampler(ExcitationGenerator const &eG_, detType const &HF,
			          Basis const &fullBasis_, Parametrization const &para_,
                int numDets_):Sampler(eG_,numDets_),para(&para_),
                lastCoeff(para_.getCoeff(HF)),fullBasis(&fullBasis_),cDet(HF){};

MetropolisSampler::MetropolisSampler(Hamiltonian const &H_, detType const &HF,
                Basis const &fullBasis_,Parametrization const &para_,
                int numDets_ ):Sampler(H_,HF,numDets_),para(&para_),
                lastCoeff(para_.getCoeff(HF)),fullBasis(&fullBasis_),cDet(HF){};

MetropolisSampler::~MetropolisSampler() {
}

//---------------------------------------------------------------------------//

void MetropolisSampler::setReference(detType const &a) {
  cDet = a;
  lastCoeff = para->getCoeff(a);
}

//---------------------------------------------------------------------------//

void MetropolisSampler::iterate(coeffType &cI, detType &dI, double &weight,
    int i){
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
  //double probRand = 0.;
  if (rng() < probRand) {
    tmp = getRandomDeterminant(*fullBasis);
    //tmp=getRandomConnection(cDet,pEx);
    // set pEx and pBack both to 1, because the back and forth prob are the same
    // and only the ratio matters.
	  pEx = 1.; 
    pBack = 1.;
  }
  else {
    tmp=Sampler::getRandomConnection(cDet,pEx);
	  pBack = Sampler::getConnectionProb(tmp,cDet);
  }
  tmpCoeff = para->getCoeff(tmp);

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
	weight =1.;
}

//---------------------------------------------------------------------------//

// create some random determinant - probably subject to change
detType MetropolisSampler::randDet() const{
	return getRandomDeterminant(*fullBasis);
}

}
