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
template <typename T=VecType>
class FullSampler: public Sampler {
public:
	FullSampler(ExcitationGenerator const &eG_, Basis const &fullBasis_,
			detType const &HF, Parametrization<T> const &para_);
	FullSampler(Hamiltonian const &H_, Basis const &fullBasis_,
			detType const &HF, Parametrization<T> const &para_);
	void iterate(coeffType &cI, detType &dI, int i) const;
	// create a dynamic polymorphic copy
	virtual FullSampler* clone() const {return new FullSampler(*this);}
	virtual ~FullSampler();
private:
	mutable int pos;
  // sampling depends on the coefficients, as they have to be given alongside the determinants
    Parametrization<T> const *para;
};

} /* namespace networkVMC */

#endif /* SRC_SAMPLERS_FULLSAMPLER_HPP_ */
