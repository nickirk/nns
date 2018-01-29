/*
 * Trainer.hpp
 *
 *  Created on: Jan 23, 2018
 *      Author: guther
 */

#ifndef SRC_TRAINER_HPP_
#define SRC_TRAINER_HPP_

#include "Nnw.hpp"
#include "CoeffType.hpp"
#include "Determinant.hpp"
#include "Sampler.hpp"

class Trainer {
public:
	Trainer(NeuralNetwork &NNW_, Sampler &msampler);
	void train(double learningRate, int method, int iteration);
	void getNextCoeff(coeffType &cI, detType &dI);
	double getE() const;
	State getState() const;
	virtual ~Trainer();
private:
	NeuralNetwork &NNW;
	Sampler &msampler;

	// containers for
	// the sampled determinants
	std::vector<detType > sampledDets;
	// their coefficients
	std::vector<coeffType > sampledCoeffs;

	// Maybe
	std::vector<std::vector<coeffType >> coupledCoeffsEpoch;
	std::vector<std::vector<detType >> coupledDetsEpoch;
};

#endif /* SRC_TRAINER_HPP_ */
