/*
 * EnergyEstimator.hpp
 * based on EnergyCF.hpp
 *  Created on: Nov 08, 2017
 *      Author: Liao
 */

#ifndef SRC_COSTFUNCTIONS_ENERGYESPREFETCHED_HPP_
#define SRC_COSTFUNCTIONS_ENERGYESPREFETCHED_HPP_

#include <vector>
#include <Eigen/Dense>
#include "../utilities/TypeDefine.hpp"
#include "CostFunction.hpp"
#include "EnergyEs.hpp"

namespace networkVMC{

class Hamiltonian;

class EnergyEsPreFetched: public CostFunction{
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
	explicit EnergyEsPreFetched(Hamiltonian const &H_):
		CostFunction(),H(H_),energy(0.0),normalizerCoeff(0.0){};

	// Hamiltonian to which the CFs energy refers
	Hamiltonian const& H;
	mutable double energy;
	mutable double normalizerCoeff;
	double evaluate(State const &input) const;

	// Has a reference member, so assignment is not a thing
	EnergyEsPreFetched& operator=(EnergyEsPreFetched const &source);
	EnergyEsPreFetched& operator=(EnergyEsPreFetched &&source);
};

}

#endif /* SRC_COSTFUNCTIONS_ENERGYESPREFETCHED_HPP_ */
