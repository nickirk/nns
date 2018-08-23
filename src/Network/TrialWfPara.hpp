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
template <typename F=std::complex<double>, typename coeffType=std::complex<double>>
class TrialWfPara: public ClonableParametrization<F, coeffType, TrialWfPara<F, coeffType> > {
  public:
    using T=Eigen::Matrix<F, Eigen::Dynamic, 1>;
	// we move the parametrizations to the TrialWfParametrization, that is, their stand-alone
	// versions are empty afterwards
	TrialWfPara(Parametrization<F, coeffType> const &basePara_, Parametrization<F, coeffType> const &trialWf_):
		basePara(basePara_.clone()),trialWf(trialWf_.clone()){};

	// return of the parameters is delegated to basePara
	T const& pars() const {return basePara->pars();}
	coeffType getBaseCoeff(detType const &det) const;
	coeffType getTrialCoeff(detType const &det) const;

	Parametrization<F, coeffType> const &base() const{return *basePara;}

	// Obtain the inner derivative dX/dPars with given dX/dC (C are coefficients)
	// delegated to basePara
	T getDeriv(detType const &det) const{return basePara->getDeriv(det);}

	// get the coefficient from base and trial para
	coeffType getCoeff(detType const &det) const;
private:
	DeepCpyUniquePtr<Parametrization<F, coeffType> > basePara;
	DeepCpyUniquePtr<Parametrization<F, coeffType> const> trialWf;
};

} /* namespace networkVMC */

#endif /* SRC_NETWORK_TRIALWFPARA_HPP_ */
