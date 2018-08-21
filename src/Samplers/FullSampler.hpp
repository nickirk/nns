/*
 * FullSampler.hpp
 *
 *  Created on: May 3, 2018
 *      Author: guther
 */

#ifndef SRC_SAMPLERS_FULLSAMPLER_HPP_
#define SRC_SAMPLERS_FULLSAMPLER_HPP_

#include "Sampler.hpp"

namespace networkVMC {

// This sampler does not sample: It loops over all determinants
// and returns them one-by-one
template <typename F=std::complex<double>, typename coeffType=std::complex<double>>
class FullSampler: public Sampler<F, coeffType> {
  public:
	FullSampler(ExcitationGenerator const &eG_, Basis const &fullBasis_,
			Parametrization<F, coeffType> const &para_);
	FullSampler(Hamiltonian const &H_, Basis const &fullBasis_,
			Parametrization<F, coeffType> const &para_);
	void iterate(coeffType &cI, detType &dI, double& weight, int i);
	// create a dynamic polymorphic copy
	virtual FullSampler* clone() const {return new FullSampler<F, coeffType>(*this);}
	virtual ~FullSampler();
private:
	int pos;
  // sampling depends on the coefficients, as they have to be given alongside the determinants
    Parametrization<F, coeffType> const *para;
    Basis const *fullBasis;
    using Sampler<F, coeffType>::numDets;
};

} /* namespace networkVMC */

#endif /* SRC_SAMPLERS_FULLSAMPLER_HPP_ */
