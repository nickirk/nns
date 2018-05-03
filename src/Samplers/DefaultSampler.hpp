/*
 * DefaultSampler.hpp
 *
 *  Created on: Apr 30, 2018
 *      Author: guther
 */

#ifndef SRC_SAMPLERS_DEFAULTSAMPLER_HPP_
#define SRC_SAMPLERS_DEFAULTSAMPLER_HPP_

#include "Sampler.hpp"

namespace networkVMC {

// This is just the Default sampler made accessible

class DefaultSampler : public Sampler{
public:
	// It works like the abstract base class, just that it can iterate
    DefaultSampler(Hamiltonian const &H_, Basis const &fullBasis_, detType const &HF, int numDets_= 100):
    	Sampler(H_,fullBasis_,HF,numDets_){};
	virtual ~DefaultSampler();
	// We use the
	virtual void iterate(coeffType &cI, detType &dI) const{
		cDet = getRandomConnection(cDet);
		dI = getDet();
		cI = coeffType();
	}
};

} /* namespace networkVMC */

#endif /* SRC_SAMPLERS_DEFAULTSAMPLER_HPP_ */
