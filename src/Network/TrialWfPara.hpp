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

// set the defaults
#include "TrialWfParaForward.hpp"

namespace networkVMC {

// this is a parametrization over a trial wavefunction (given as a parametrization, too)
template <typename F, typename coeffType>
class TrialWfPara: public ClonableParametrization<F, coeffType, TrialWfPara<F, coeffType> > {
  public:
    using T=Eigen::Matrix<F, Eigen::Dynamic, 1>;
	// we move the parametrizations to the TrialWfParametrization, that is, their stand-alone
	// versions are empty afterwards
	TrialWfPara(Parametrization<F, coeffType> const &basePara_, Parametrization<F, coeffType> const &trialWf_):
		// copy the parametrizations to the stored ones
		basePara(basePara_.clone()),trialWf(trialWf_.clone()){};

	// moving version of construction
	template<typename P1, typename P2>
	static TrialWfPara<F, coeffType> toTrialWfPara(P1&& basePara_, P2&& trialWf_){
		return TrialWfPara(basePara_.moveClone(), trialWf_.moveClone());
	}

	// return of the parameters is delegated to basePara
	T const& pars() const {return basePara->pars();}
	coeffType getBaseCoeff(detType const &det) const;

	// additional functionality of the trialwfpara - return the coeffs of the trial wf
	coeffType getTrialCoeff(detType const &det) const;

	Parametrization<F, coeffType> const &base() const{return *basePara;}

	// Obtain the inner derivative dX/dPars with given dX/dC (C are coefficients)
	// delegated to basePara
	T getDeriv(detType const &det) const{return basePara->getDeriv(det);}
	// this is NOT the same as this->Parametrization<F, coeffType>::getMarkovDeriv !!!
	T getMarkovDeriv(detType const &det) const{return basePara->getMarkovDeriv(det);}

	// get the coefficient from base and trial para
	coeffType getCoeff(detType const &det) const;
private:
	DeepCpyUniquePtr<Parametrization<F, coeffType> > basePara;
	DeepCpyUniquePtr<Parametrization<F, coeffType> const> trialWf;

	// auxiliary ctor for ownership-transferring construction
	// this is too dangerous to be used outside the class
	TrialWfPara(Parametrization<F, coeffType> *base, Parametrization<F, coeffType> *trial):
		basePara(base),trialWf(trial){};  // initialize with given pointers
};

} /* namespace networkVMC */

#endif /* SRC_NETWORK_TRIALWFPARA_HPP_ */
