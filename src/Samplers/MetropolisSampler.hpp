/*
 * MarkovSampler.hpp
 *
 *  Created on: Jan 23, 2018
 *      Author: guther
 */

#ifndef SRC_SAMPLERS_METROPOLISSAMPLER_HPP_
#define SRC_SAMPLERS_METROPOLISSAMPLER_HPP_

#include "../HilbertSpace/Determinant.hpp"
#include "../Network/Nnw.hpp"
#include "Sampler.hpp"

// Class for sampling determinants using the metropolis algorithm

class MetropolisSampler: public Sampler {
public:
	MetropolisSampler(Hamiltonian const &H_, Basis const &fullBasis_, detType const &HF, NeuralNetwork const &NNW_):
			Sampler(H_,fullBasis_,HF),NNW(NNW_),lastCoeff(NNW_.getCoeff(cDet)){};
	virtual ~MetropolisSampler();
	//Do a markov step
	virtual void iterate(coeffType &cI, detType &dI) const;
private:
	NeuralNetwork const &NNW;
	mutable coeffType lastCoeff;
};

#endif /* SRC_SAMPLERS_METROPOLISSAMPLER_HPP_ */
