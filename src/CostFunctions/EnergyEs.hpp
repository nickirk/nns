/*
 * EnergyEs.hpp
 *
 *  Created on: Jun 26, 2018
 *      Author: guther
 */

#ifndef SRC_COSTFUNCTIONS_ENERGYES_HPP_
#define SRC_COSTFUNCTIONS_ENERGYES_HPP_

#include "CostFunction.hpp"
#include "EnergyCFBaseClass.hpp"
#include "../utilities/DeepCpyUniquePtr.hpp"

namespace networkVMC {

class Hamiltonian;

class EnergyEs: public CostFunction {
public:
	// uses the required number of connections and the Hamiltonian
	EnergyEs(Hamiltonian const &H_, int numCons_=20);
	virtual ~EnergyEs();
	virtual void setUpCF(SamplerType const &sT);

// The operations are actually performed by another EnergyCF (worker)
	nablaType nabla(State const &input) const {return worker->nabla(input);}
	double calc(State const &input) const {return worker->calc(input);};

	double getNormalizer() const {return worker->getNormalizer();}

	// Allow for polymorphic copy
	virtual EnergyEs* clone() const {return new EnergyEs(*this);}

	// the energy estimators do need connections
	virtual int connectionsRequired() const {return numCons;}
private:
	Hamiltonian const& H;
	// this is not a stand-alone CF, the work is done by another,
	// owned CF
	DeepCpyUniquePtr<EnergyCFBaseClass> worker;
	// number of connected states to be taken into account
	int numCons;
};

} /* namespace networkVMC */

#endif /* SRC_COSTFUNCTIONS_ENERGYES_HPP_ */
