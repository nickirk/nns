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

namespace networkVMC{

class Hamiltonian;

class EnergyEsMarkov: public CostFunction{
public:
	explicit EnergyEsMarkov(Hamiltonian const &H_):CostFunction(),H(H_),energy(0.0),normalizerCoeff(0.0){};
	nablaType nabla(State const &input) const;
	double calc(State const &input) const {return energy;};
    double getNormalizer(){return normalizerCoeff;};
private:
	Hamiltonian const& H;
	mutable double energy;
	mutable double normalizerCoeff;
	double evaluate(State const &input) const;

	// Has a reference member, so assignment is not a thing
	EnergyEsMarkov& operator=(EnergyEsMarkov const &source);
	EnergyEsMarkov& operator=(EnergyEsMarkov &&source);
};

}

#endif /* SRC_COSTFUNCTIONS_ENERGYESMARKOV_HPP_ */
