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
#include "Eigen/Dense"
#include "TypeDefine.hpp"
// This cost function tries to minimize the energy expectation value

class EnergyCF: public CostFunction{
public:
	// Here, we need to supply a Hamiltonian
	explicit EnergyCF(Hamiltonian const &H_):CostFunction(),H(H_),energy(0.0),normalizerCoeff(0.0){};
// implementation of the function itself and its derivative
	virtual std::vector<Eigen::VectorXd> nabla(std::vector<State> const &input) const;
	virtual double calc(std::vector<State> const &input) const {return energy;}
// auxiliary function to compute the denominator of the expectation value
    double getNormalizer(){return normalizerCoeff;};
private:
	Hamiltonian const& H;
// treat energy and normalizer as mutable, as they are only auxiliary cache variables
	mutable double energy;
	mutable double normalizerCoeff;
	double evaluate(std::vector<State> const &input) const;
};



#endif /* SRC_ENERGYCF_HPP_ */
