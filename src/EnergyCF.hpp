/*
 * EnergyCF.hpp
 *
 *  Created on: Oct 25, 2017
 *      Author: guther
 */

#ifndef SRC_ENERGYCF_HPP_
#define SRC_ENERGYCF_HPP_

#include <vector>
#include "CostFunction.hpp"
#include "State.hpp"
#include "Hamiltonian.hpp"
#include "CoeffType.hpp"

class EnergyCF: public CostFunction{
public:
	explicit EnergyCF(Hamiltonian const &H_):CostFunction(),H(H_),energy(0.0),normalizerCoeff(0.0){};
	friend class Evaluator;
private:
	Hamiltonian const& H;
	mutable double energy;
	mutable double normalizerCoeff;
	std::vector<coeffType > nabla(State const &input) const;
	void calc(State const &input) const;
	double getValue() const{return energy;}
};



#endif /* SRC_ENERGYCF_HPP_ */
