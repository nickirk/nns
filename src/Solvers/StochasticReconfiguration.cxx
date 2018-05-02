/*
 * StochasticReconfiguration.cxx
 *
 *  Created on: May 2, 2018
 *      Author: guther
 */

#include "StochasticReconfiguration.hpp"

namespace networkVMC {

StochasticReconfiguration::~StochasticReconfiguration() {
}

//---------------------------------------------------------------------------------------------------//

// Conjugate descent (I ported this from Nnw.cxx and will not touch it)
void StochasticReconfiguration::update(VecType &w, VecType const &force,
		  State const &input) const{

	// first, get the input vector's coefficients
	std::size_t numDets = input.size();
	VecType ci(numDets);
	for(std::size_t i = 0; i<numDets; ++ i){
		ci[i] = input.coeff(i);
	}

	// and the second order derivative
	auto dcdw = par.calcdCdwSR(input);

	Eigen::MatrixXd s, okokp;
	Eigen::VectorXcd ok;
	double normalizer = ci.norm();
	ok = dcdw*ci.conjugate()/normalizer;
	okokp = (dcdw*dcdw.adjoint()/normalizer).real();
	s = okokp - (ok*ok.adjoint()).real();
        double lambda = std::max(100*std::pow(0.999,iteration), 1e-2);
        s+=s.diagonal().asDiagonal()*lambda;
	w-=gamma*s.inverse()*force;
	// increase the iteration counter
	iteration += 1;
}

} /* namespace networkVMC */
