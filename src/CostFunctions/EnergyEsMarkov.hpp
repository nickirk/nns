/*
 * EnergyEsMarkov.hpp
 * based on EnergyCF.hpp
 *  Created on: Jan 29, 2018
 *      Author: Liao
 */

#ifndef SRC_COSTFUNCTIONS_ENERGYESMARKOV_HPP_
#define SRC_COSTFUNCTIONS_ENERGYESMARKOV_HPP_

#include <vector>
#include <Eigen/Dense>
#include "../utilities/TypeDefine.hpp"
#include "CostFunction.hpp"
#include "../Hamiltonian/Hamiltonian.hpp"
#include "../utilities/State.hpp"

class EnergyEsMarkov: public CostFunction{
public:
	explicit EnergyEsMarkov(Hamiltonian const &H_):CostFunction(),H(H_),energy(0.0),normalizerCoeff(0.0){};
	std::vector<Eigen::VectorXd> nabla(std::vector<State > const &input) const;
	double calc(std::vector<State > const &input) const {return energy;};
        double getNormalizer(){return normalizerCoeff;};
private:
	Hamiltonian const& H;
	mutable double energy;
	mutable double normalizerCoeff;
	double evaluate(std::vector<State > const &input) const;
};



#endif /* SRC_COSTFUNCTIONS_ENERGYESMARKOV_HPP_ */
