/*
 * TrialWfParametrization.h
 *
 *  Created on: Aug 16, 2018
 *      Author: guther
 */

#ifndef SRC_NETWORK_TRIALWFPARA_HPP_
#define SRC_NETWORK_TRIALWFPARA_HPP_

#include "Parametrization.hpp"
#include "../utilities/DeepCpyUniquePtr.hpp"

namespace networkVMC {

// this is a parametrization over a trial wavefunction (given as a parametrization, too)
template<typename T=VecType>
class TrialWfPara: public ClonableParametrization<T,TrialWfPara<T> > {
public:
	// we move the parametrizations to the TrialWfParametrization, that is, their stand-alone
	// versions are empty afterwards
	TrialWfPara(Parametrization<T> const &basePara_, Parametrization<T> const &trialWf_):
		basePara(basePara_.clone()),trialWf(trialWf_.clone()){};
	virtual ~TrialWfPara() = default;

	// return of the parameters is delegated to basePara
	virtual T const& pars() const {return basePara->pars();}

	// get the coefficient
	virtual coeffType getCoeff(detType const &det) const;
	// Obtain the inner derivative dX/dPars with given dX/dC (C are coefficients)
	// delegated to basePara
	virtual T calcNablaPars(State const &input, nablaType const &outerDerivative);
private:
	DeepCpyUniquePtr<Parametrization<T> > basePara;
	DeepCpyUniquePtr<Parametrization<T> const> trialWf;
};

} /* namespace networkVMC */

#endif /* SRC_NETWORK_TRIALWFPARA_HPP_ */
