/*
 * EnergyEstimator.hpp
 * based on EnergyCF.hpp
 *  Created on: Nov 08, 2017
 *      Author: Liao
 */

#ifndef SRC_COSTFUNCTIONS_ENERGYESTIMATOR_HPP_
#define SRC_COSTFUNCTIONS_ENERGYESTIMATOR_HPP_

#include <vector>
#include <Eigen/Dense>
#include "../utilities/TypeDefine.hpp"
#include "CostFunction.hpp"

namespace networkVMC{

class Hamiltonian;

class EnergyEstimator: public CostFunction{
public:
	explicit EnergyEstimator(Hamiltonian const &H_):CostFunction(),H(H_),energy(0.0),normalizerCoeff(0.0){};
	std::vector<Eigen::VectorXd> nabla(State const &input) const;
	double calc(State const &input) const {return energy;};
    double getNormalizer(){return normalizerCoeff;};
private:
	Hamiltonian const& H;
	mutable double energy;
	mutable double normalizerCoeff;
	double evaluate(State const &input) const;

	// Has a reference member, so assignment is not a thing
	EnergyEstimator& operator=(EnergyEstimator const &source);
	EnergyEstimator& operator=(EnergyEstimator &&source);
};

}

#endif /* SRC_COSTFUNCTIONS_ENERGYESTIMATOR_HPP_ */
