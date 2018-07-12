/*
 * StochasticReconfiguration.cxx
 *
 *  Created on: May 2, 2018
 *      Author: guther
 */

#include "StochasticReconfiguration.hpp"

namespace networkVMC {

template <typename T>
StochasticReconfiguration<T>::~StochasticReconfiguration() {
}

//---------------------------------------------------------------------------------------------------//

// Conjugate descent (I ported this from Nnw.cxx and will not touch it)
template <typename T>
void StochasticReconfiguration<T>::update(T &w, T const &force,
		  State const &input){

	// first, get the input vector's coefficients
	std::size_t numDets = input.size();
	Eigen::VectorXcd ci(numDets);
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
    double lambda = std::max(10*std::pow(0.999,iteration), 1.);
    s+=s.diagonal().asDiagonal()*lambda;
	w-=Solver<T>::learningRate*s.inverse()*force;
	// increase the iteration counter
	iteration += 1;
}
template class StochasticReconfiguration<VecType>;
template class StochasticReconfiguration<VecCType>;
} /* namespace networkVMC */
