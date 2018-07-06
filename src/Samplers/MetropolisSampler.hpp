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
	MetropolisSampler(ExcitationGenerator const &eG_, Basis const &fullBasis_, detType const &HF, Parametrization<T> const &para_):
			Sampler(eG_,fullBasis_,HF),para(&para_),lastCoeff(para_.getCoeff(cDet)){};
	MetropolisSampler(Hamiltonian const &H_, Basis const &fullBasis_, detType const &HF, Parametrization<T> const &para_):
			Sampler(H_,fullBasis_,HF),para(&para_),lastCoeff(para_.getCoeff(cDet)){};
	virtual ~MetropolisSampler();
	// create a dynamic polymorphic copy
	virtual MetropolisSampler* clone() const {return new MetropolisSampler(*this);}
	//Do a markov step
	virtual void iterate(coeffType &cI, detType &dI, double &weight, int i) const;
	// this is a markov-type sampler
	SamplerType type() const {return Markov;}
private:
  // sampling depends on the coefficients, as they have to be given alongside the determinants
    Parametrization<T> const *para;
	mutable coeffType lastCoeff;
};

}
#endif /* SRC_SAMPLERS_METROPOLISSAMPLER_HPP_ */
