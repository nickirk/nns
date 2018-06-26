/*
 * EnergyEs.hpp
 *
 *  Created on: Jun 26, 2018
 *      Author: guther
 */

#ifndef SRC_COSTFUNCTIONS_ENERGYES_HPP_
#define SRC_COSTFUNCTIONS_ENERGYES_HPP_

#include "CostFunction.hpp"

namespace networkVMC {

class Hamiltonian;

class EnergyEs: public CostFunction {
public:
	EnergyEs(Hamiltonian const &H_);
	virtual ~EnergyEs();
	virtual std::unique_ptr<CostFunction> setUpCF(SamplerType const &sT) const;

// Even though this is never gonna be used, it has to be defined
// TODO route these definitions via the correct EnergyEs
	nablaType nabla(State const &input) const {return nablaType();}
	double calc(State const &input) const {return 0.0;};

	// Allow for polymorphic copy
	virtual CostFunction* clone() const {return new EnergyEs(*this);}
private:
	Hamiltonian const& H;
};

} /* namespace networkVMC */

#endif /* SRC_COSTFUNCTIONS_ENERGYES_HPP_ */
