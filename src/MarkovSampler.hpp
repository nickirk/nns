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
	MarkovSampler(Hamiltonian const &H_, Basis const &fullBasis_, int numDets_, detType const &HF):
			Sampler(H_,fullBasis_,numDets_,HF){};
	virtual ~MarkovSampler();
	//Do a markov step
	void iterate(coeffType &cI, detType &dI, NeuralNetwork const &NNW) const;
};

#endif /* SRC_MARKOVSAMPLER_HPP_ */
