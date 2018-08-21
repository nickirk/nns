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
template <typename F=std::complex<double>, typename coeffType=std::complex<double>>
class DefaultSampler : public Sampler{
  public:
	// It works like the abstract base class, just that it can iterate
    DefaultSampler(ExcitationGenerator const &eG_, Basis const &fullBasis_, detType const &HF,
    		Parametrization<F, coeffType> const &para_, int numDets_= 100):
    	Sampler(eG_,HF,numDets_),para(&para_),fullBasis(&fullBasis_){};

	// It works like the abstract base class, just that it can iterate
    DefaultSampler(Hamiltonian const &H_, Basis const &fullBasis_, detType const &HF,
    		Parametrization<F, coeffType> const &para_, int numDets_= 100):
    	Sampler(H_,HF,numDets_),para(&para_),fullBasis(&fullBasis){};

	virtual ~DefaultSampler();
	// create a dynamic polymorphic copy
	virtual DefaultSampler* clone() const {return new DefaultSampler(*this);}
	// We use the default functionality for iterate()
	// this class just makes it accessible for testing purposes
	virtual void iterate(F &cI, detType &dI, double& weight, int i){
        //TODO needs to fill up the weight
		double p;
		cDet = getRandomConnection(cDet,p);
		dI = getDet();
		cI = para->getCoeff(cDet);
	}
private:
	  // sampling depends on the coefficients, as they have to be given alongside the determinants
	  Parametrization<F, coeffType> const *para;
	  // and the corresponding basis including the information on the number of electrons
	  // with a given spin
	  Basis const *fullBasis;
};

} /* namespace networkVMC */

#endif /* SRC_SAMPLERS_DEFAULTSAMPLER_HPP_ */
