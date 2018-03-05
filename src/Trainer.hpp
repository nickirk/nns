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
#include "State.hpp"

class Trainer {
public:
	Trainer(Hamiltonian const &modelHam_,NeuralNetwork &NNW_, Sampler &msampler);
	void train(double learningRate, int method, int iteration);
	void getNextCoeff(coeffType &cI, detType &dI);
	double getE() const;
	std::vector<State> getState() const;
	virtual ~Trainer();
private:
	NeuralNetwork &NNW;
	Hamiltonian const &modelHam;
	Sampler &msampler;
	std::vector<State > inputState;

};

#endif /* SRC_TRAINER_HPP_ */
