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
#include "EnergyCFBaseClass.hpp"

namespace networkVMC{

class Hamiltonian;

class EnergyEsMarkov: public EnergyCFBaseClass{
public:
	// Do it like this: We cannot generate EnergyEsMarkov
	// directly, this has to be done via EnergyEs
	friend EnergyEs;

	virtual double calc(State const &input) const {return energy;};
	nablaType nabla(State const &input) const;

	// Allow for polymorphic copy
	virtual EnergyEsMarkov* clone() const {return new EnergyEsMarkov(*this);}
private:
    // Make sure this is not manually constructed, but only via
    // EnergyEs. This way, we cannot attribute the wrong CF to a sampler
	explicit EnergyEsMarkov(Hamiltonian const &H_):
		EnergyCFBaseClass(H_){};

	double evaluate(State const &input) const;
};

}

#endif /* SRC_COSTFUNCTIONS_ENERGYESMARKOV_HPP_ */
