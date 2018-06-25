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
#include "EnergyEs.hpp"

namespace networkVMC{

class Hamiltonian;

class EnergyEsMarkov: public CostFunction{
public:
	// Do it like this: We cannot generate EnergyEsMarkov
	// directly, this has to be done via EnergyEs
	friend EnergyEs;

	nablaType nabla(State const &input) const;
	double calc(State const &input) const {return energy;};
    double getNormalizer(){return normalizerCoeff;};
private:
    // Make sure this is not manually constructed, but only via
    // EnergyEs. This way, we cannot attribute the wrong CF to a sampler
	explicit EnergyEsMarkov(Hamiltonian const &H_):
		CostFunction(),H(H_),energy(0.0),normalizerCoeff(0.0){};

	// Hamiltonian, to which the energy referes
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
