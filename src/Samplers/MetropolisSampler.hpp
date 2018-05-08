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

class MetropolisSampler: public Sampler {
public:
	MetropolisSampler(Hamiltonian const &H_, Basis const &fullBasis_, detType const &HF, Parametrization const &para_):
			Sampler(H_,fullBasis_,HF,para_),lastCoeff(para_.getCoeff(cDet)){};
	virtual ~MetropolisSampler();
	//Do a markov step
	virtual void iterate(coeffType &cI, detType &dI) const;
private:
	mutable coeffType lastCoeff;
};

}
#endif /* SRC_SAMPLERS_METROPOLISSAMPLER_HPP_ */
