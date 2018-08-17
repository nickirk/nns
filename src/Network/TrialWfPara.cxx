/*
 * TrialWfParametrization.cpp
 *
 *  Created on: Aug 16, 2018
 *      Author: guther
 */

#include "TrialWfPara.hpp"
#include "../utilities/State.hpp"

namespace networkVMC {

template<typename T>
coeffType TrialWfPara<T>::getCoeff(detType const &det) const{
	// the coefficient is the product of trial wf and parametrized coeff
	return basePara->getCoeff(det)*trialWf->getCoeff(det);
}

template<typename T>
T TrialWfPara<T>::calcNablaPars(State const &input, nablaType const &outerDerivative){
	nablaType scaledOuterDerivative = outerDerivative;
	// scale the outerDerivative with the trial WF

	// number of dets in input
	size_t numDets = input.size();
	// total number of entries in input, including coupled dets
	size_t totalDets = input.totalSize();

	// first, the derivative w.r. to sampled coeffs
	for(size_t i = 0; i < numDets; ++i){
		scaledOuterDerivative[i] /= trialWf->getCoeff(input.det(i));
	}

	// then w.r. to the connected coeffs
	int pos=0;
	std::vector<coeffType> cDets;
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

// instantiate template classes
template class TrialWfPara<VecType>;
template class TrialWfPara<VecCType>;
} /* namespace networkVMC */
