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
/*
template <typename F, typename coeffType>
TrialWfPara<F, coeffType>::T TrialWfPara<F, coeffType>::calcNablaPars(State<coeffType> const &input, T const &outerDerivative){
	T scaledOuterDerivative = outerDerivative;
	// scale the outerDerivative with the trial WF

	// number of dets in input
	size_t numDets = input.size();
	// total number of entries in input, including coupled dets
	size_t totalDets = input.totalSize();

	// the outer derivative is the derivative with respect to all parameters of the State input
	// so those two numbers better be the same
	if(outerDerivative.size() != totalDets) throw SizeMismatchError(outerDerivative.size(),totalDets);

	// first, the derivative w.r. to sampled coeffs
	for(size_t i = 0; i < numDets; ++i){
		scaledOuterDerivative[i] /= trialWf->getCoeff(input.det(i));
	}

	// then w.r. to the connected coeffs
	int pos=0;
	std::vector<detType> cDets;
	for(size_t i = 0; i < numDets; ++i){
		cDets = input.coupledDets(i);
		for(size_t j = 0; j < cDets.size(); ++j){
			pos = input.locate(i);
			scaledOuterDerivative[j+numDets+pos] /= trialWf->getCoeff(cDets[j]);
		}
	}

	// then, add the inner derivative
	return basePara->calcNablaPars(input, scaledOuterDerivative);
}
*/
// instantiate template classes
template class TrialWfPara<double, double>;
;
template class TrialWfPara<std::complex<double>, std::complex<double>>;
} /* namespace networkVMC */
