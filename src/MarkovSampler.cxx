/*
 * MarkovSampler.cxx
 *
 *  Created on: Jan 23, 2018
 *      Author: guther
 */

#include "MarkovSampler.hpp"
#include "Hamiltonian.hpp"
#include <random>
#include <cmath>

MarkovSampler::~MarkovSampler() {
}

void MarkovSampler::iterate(coeffType &cI, detType &dI) const{
	int const maxIters{1000};
	double p{0.0};

	// This is tricky, we need to cache the network state
	// but we don't know, which one at this point
	// NNW.cacheNetworkState();
	// set up the rng
	std::random_device rd;     // only used once to initialise (seed) engine
	std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
	std::uniform_real_distribution<double> uni;		// Uniform distribution from 0.0 to 1.0
	// Iterate until acceptance (or maximum number of iters)
	//for(int i=0;i<maxIters;++i){
	// First, get a random coupled determinant (from cDet)
        detType tmp{getRandomCoupledState(cDet, p)};
	// Not a coupled, but a total random determinant.
	//std::vector<int > spinConfig=fullBasis.getSpinConfig();
	//detType tmp{getRandomDeterminant(spinConfig)};
	// And its coefficient
	coeffType tmpCoeff{NNW.getCoeff(tmp)};
	if(uni(rng) - std::norm(tmpCoeff)/std::norm(lastCoeff)<-1e-8){
		// With probability |cJ|^2/|cI|^2, accept the move
		cDet = tmp;
		lastCoeff = tmpCoeff;
		// We stored the state before, but the wrong one. Correct this
		NNW.cacheNetworkState();
	}
        else {
         NNW.repCachedNetworkState(); 
        }

	// and set the output
	cI = lastCoeff;
	dI = cDet;
}
