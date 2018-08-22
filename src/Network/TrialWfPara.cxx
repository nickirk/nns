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

template<typename T>
coeffType TrialWfPara<T>::getCoeff(detType const &det) const{
	// the coefficient is the product of trial wf and parametrized coeff
	return basePara->getCoeff(det)*trialWf->getCoeff(det);
}

//---------------------------------------------------------------------------//

template<typename T>
T TrialWfPara<T>::calcNablaParsConnected(State const &input, nablaType const &outerDerivative){
	auto scaledOuterDerivative = scaleDerivative(input,outerDerivative);

	// then, add the inner derivative
	return basePara->calcNablaParsConnected(input, scaledOuterDerivative);
}

//---------------------------------------------------------------------------//

template<typename T>
T TrialWfPara<T>::calcNablaParsMarkovConnected(State const &input, nablaType const &outerDerivative){
	auto scaledOuterDerivative = scaleDerivative(input,outerDerivative);

	// then, add the inner derivative
	return basePara->calcNablaParsMarkovConnected(input, scaledOuterDerivative);
}

//---------------------------------------------------------------------------//

template<typename T>
nablaType TrialWfPara<T>::scaleDerivative(State const &input, nablaType const &outerDerivative) const{
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
		outerDerivative[i] *= trialWf->getCoeff(input.det(i));
	}

	// then w.r. to the connected coeffs
	int pos=0;
	std::vector<detType> cDets;
	for(size_t i = 0; i < numDets; ++i){
		cDets = input.coupledDets(i);
		for(size_t j = 0; j < cDets.size(); ++j){
			pos = input.locate(i);
			outerDerivative[j+numDets+pos] *= trialWf->getCoeff(cDets[j]);
		}
	}

	return outerDerivative;
}

// instantiate template classes
template class TrialWfPara<VecType>;
template class TrialWfPara<VecCType>;
} /* namespace networkVMC */
