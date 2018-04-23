/*
 * MarkovSampler.hpp
 *
 *  Created on: Jan 23, 2018
 *      Author: guther
 */

#ifndef SRC_METROPOLISSAMPLER_HPP_
#define SRC_METROPOLISSAMPLER_HPP_

#include "Sampler.hpp"
#include "Nnw.hpp"
#include "Determinant.hpp"

// Class for sampling determinants using the metropolis algorithm

class MetropolisSampler: public Sampler {
public:
	MetropolisSampler(Hamiltonian const &H_, Basis const &fullBasis_, int numDets_, detType const &HF, NeuralNetwork const &NNW_):
			Sampler(H_,fullBasis_,numDets_,HF),NNW(NNW_),lastCoeff(NNW_.getCoeff(cDet)){};
	virtual ~MetropolisSampler();
	//Do a markov step
	void iterate(coeffType &cI, detType &dI) const;
private:
	NeuralNetwork const &NNW;
	mutable coeffType lastCoeff;
};

#endif /* SRC_METROPOLISSAMPLER_HPP_ */
