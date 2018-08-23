/*
 * TrialWfParametrization.cpp
 *
 *  Created on: Aug 16, 2018
 *      Author: guther
 */

#include "TrialWfPara.hpp"
#include "../utilities/State.hpp"
#include "../utilities/Errors.hpp"

namespace networkVMC {

template <typename F, typename coeffType>
coeffType TrialWfPara<F, coeffType>::getCoeff(detType const &det) const{
	// the coefficient is the product of trial wf and parametrized coeff
	return basePara->getCoeff(det)*trialWf->getCoeff(det);
}

template <typename F, typename coeffType>
coeffType TrialWfPara<F, coeffType>::getBaseCoeff(detType const &det) const{
	return basePara->getCoeff(det);
}

template <typename F, typename coeffType>
coeffType TrialWfPara<F, coeffType>::getTrialCoeff(detType const &det) const{
	return trialWf->getCoeff(det);
}
// instantiate template classes
template class TrialWfPara<double, double>;
template class TrialWfPara<std::complex<double>, std::complex<double>>;
} /* namespace networkVMC */
