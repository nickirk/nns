/*
 * EnergyEstimator.hpp
 * based on EnergyCF.hpp
 *  Created on: Nov 08, 2017
 *      Author: Liao
 */

#ifndef SRC_ENERGYESTIMATOR_HPP_
#define SRC_ENERGYESTIMATOR_HPP_

#include <vector>
#include <Eigen/Dense>
#include "CostFunction.hpp"
#include "State.hpp"
#include "Hamiltonian.hpp"
#include "CoeffType.hpp"

class EnergyEstimator: public CostFunction{
public:
	explicit EnergyEstimator(Hamiltonian const &H_):CostFunction(),H(H_),energy(0.0),normalizerCoeff(0.0){};
	std::vector<Eigen::VectorXd> nabla(std::vector<State> const &input) const;
	double calc(std::vector<State> const &input) const {return energy;};
        double getNormalizer(){return normalizerCoeff;};
private:
	Hamiltonian const& H;
	mutable double energy;
	mutable double normalizerCoeff;
	double evaluate(std::vector<State> const &input) const;
};



#endif /* SRC_ENERGYESTIMATOR_HPP_ */
