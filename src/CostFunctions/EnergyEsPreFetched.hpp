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
#include "EnergyCFBaseClass.hpp"

namespace networkVMC{

class Hamiltonian;

class EnergyEsPreFetched: public EnergyCFBaseClass{
public:
	// Do it like this: We cannot generate EnergyEsPreFetched
	// directly, this has to be done via EnergyEs
	friend EnergyEs;

	nablaType nabla(State const &input) const;

	// Allow for polymorphic copy
	virtual EnergyEsPreFetched* clone() const {return new EnergyEsPreFetched(*this);}

	// For sake of completeness, we specify that this requires connections
	virtual int connectionsRequired() const {return numCons;}
private:
    // Make sure this is not manually constructed, but only via
    // EnergyEs. This way, we cannot attribute the wrong CF to a sampler
	explicit EnergyEsPreFetched(Hamiltonian const &H_, int numCons_):
		EnergyCFBaseClass(H_), numCons(numCons_){};

	coeffType evaluate(State const &input) const;

	int numCons;
};

}

#endif /* SRC_COSTFUNCTIONS_ENERGYESPREFETCHED_HPP_ */
