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
template <typename T=VecType>
class DefaultSampler : public Sampler{
public:
	// It works like the abstract base class, just that it can iterate
    DefaultSampler(Hamiltonian const &H_, Basis const &fullBasis_, detType const &HF,
    		Parametrization<T> const &para_, int numDets_= 100):
    	Sampler(H_,fullBasis_,HF,numDets_),para(&para_){};
	virtual ~DefaultSampler();
	// create a dynamic polymorphic copy
	virtual DefaultSampler* clone() const {return new DefaultSampler(*this);}
	// We use the default functionality for iterate()
	// this class just makes it accessible for testing purposes
	virtual void iterate(coeffType &cI, detType &dI, int i) const{
		double p;
		cDet = getRandomConnection(cDet,p);
		dI = getDet();
		cI = para->getCoeff(cDet);
	}
private:
	  // sampling depends on the coefficients, as they have to be given alongside the determinants
	    Parametrization<T> const *para;

};

} /* namespace networkVMC */

#endif /* SRC_SAMPLERS_DEFAULTSAMPLER_HPP_ */
