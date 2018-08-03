/*
 * MarkovSampler.hpp
 *
 *  Created on: Jan 23, 2018
 *      Author: guther
 */

#ifndef SRC_SAMPLERS_METROPOLISSAMPLER_HPP_
#define SRC_SAMPLERS_METROPOLISSAMPLER_HPP_

#include "Sampler.hpp"
#include "../Network/Parametrization.hpp"

// Class for sampling determinants using the metropolis algorithm

namespace networkVMC{
template <typename T=VecType>
class MetropolisSampler: public Sampler {
public:
	MetropolisSampler(ExcitationGenerator const &eG_, detType const &HF,
			          Basis const &fullBasis_, Parametrization<T> const &para_, 
                int numDets_ = 100):Sampler(eG_,HF,numDets_),para(&para_),
                lastCoeff(para_.getCoeff(cDet)),fullBasis(&fullBasis_){};
	MetropolisSampler(Hamiltonian const &H_, detType const &HF, 
                Basis const &fullBasis_,Parametrization<T> const &para_, 
                int numDets_ = 100):Sampler(H_,HF,numDets_),para(&para_),
                fullBasis(&fullBasis_),lastCoeff(para_.getCoeff(cDet)){};
	virtual ~MetropolisSampler();
	// create a dynamic polymorphic copy
	virtual MetropolisSampler* clone() const {return new MetropolisSampler(*this);}
	//Do a markov step
	virtual void iterate(coeffType &cI, detType &dI, double &weight, int i);
	// this is a markov-type sampler
	SamplerType type() const {return Markov;}

	void resetSpecs(){lastCoeff = para->getCoeff(cDet);}
private:
  // sampling depends on the coefficients, as they have to be given alongside the determinants
  Parametrization<T> const *para;
	coeffType lastCoeff;
  Basis const *fullBasis;
};

}
#endif /* SRC_SAMPLERS_METROPOLISSAMPLER_HPP_ */
