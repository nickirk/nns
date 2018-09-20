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

coeffType TrialWfPara::getCoeff(detType const &det) const{
	// the coefficient is the product of trial wf and parametrized coeff
	return basePara->getCoeff(det)*trialWf->getCoeff(det);
}

//---------------------------------------------------------------------------//

coeffType TrialWfPara::getBaseCoeff(detType const &det) const{
	return basePara->getCoeff(det);
}

//---------------------------------------------------------------------------//

coeffType TrialWfPara::getTrialCoeff(detType const &det) const{
	return trialWf->getCoeff(det);
}

//---------------------------------------------------------------------------//
} /* namespace networkVMC */
