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
	double p{0.0};
	// First, get a random coupled determinant (from cDet)
	detType tmp{getRandomCoupledState(cDet,p)};
	// And its coefficient
	coeffType tmpCoeff{NNW.getCoeff(tmp)};
	// set up the rng
	std::random_device rd;     // only used once to initialise (seed) engine
	std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
	std::uniform_real_distribution<double> uni;		// Uniform distribution from 0.0 to 1.0
	if(uni(rng) < std::fabs(tmpCoeff/lastCoeff)){
		// With probability cJ/cI, accept the move
		cDet = tmp;
		lastCoeff = tmpCoeff;
	}
	// and set the output
	cI = lastCoeff;
	dI = cDet;
}
