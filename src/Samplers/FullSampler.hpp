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

class Parametrization;
// This sampler does not sample: It loops over all determinants
// and returns them one-by-one
class FullSampler: public Sampler {
  public:
	FullSampler(ExcitationGenerator const &eG_, Basis const &fullBasis_,
			Parametrization const &para_);
	FullSampler(Hamiltonian const &H_, Basis const &fullBasis_,
			Parametrization const &para_);
	void iterate(coeffType &cI, detType &dI, double& weight, int i);
	// create a dynamic polymorphic copy
	virtual FullSampler* clone() const {return new FullSampler(*this);}
	virtual ~FullSampler();
private:
	int pos;
  // sampling depends on the coefficients, as they have to be given alongside the determinants
    Parametrization const *para;
    Basis const *fullBasis;
    using Sampler::numDets;
};

} /* namespace networkVMC */

#endif /* SRC_SAMPLERS_FULLSAMPLER_HPP_ */
