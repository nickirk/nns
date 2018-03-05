/*
 * MarkovSampler.hpp
 *
 *  Created on: Jan 23, 2018
 *      Author: guther
 */

#ifndef SRC_MARKOVSAMPLER_HPP_
#define SRC_MARKOVSAMPLER_HPP_

#include "Sampler.hpp"
#include "Nnw.hpp"
#include "Determinant.hpp"

class MarkovSampler: public Sampler {
public:
	MarkovSampler(Hamiltonian const &H_, Basis const &fullBasis_, int numDets_, detType const &HF, NeuralNetwork const &NNW_):
			Sampler(H_,fullBasis_,numDets_,HF),NNW(NNW_),lastCoeff(NNW_.getCoeff(cDet)){};
	virtual ~MarkovSampler();
	//Do a markov step
	void iterate(coeffType &cI, detType &dI) const;
private:
	NeuralNetwork const &NNW;
	mutable coeffType lastCoeff;
};

#endif /* SRC_MARKOVSAMPLER_HPP_ */
