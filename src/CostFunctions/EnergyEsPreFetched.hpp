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
#include "EnergyEsForward.hpp"
#include "EnergyCFBaseClass.hpp"

namespace networkVMC{

class Hamiltonian;
template <typename F=std::complex<double>, typename coeffType=std::complex<double>>
class EnergyEsPreFetched: public EnergyCFBaseClass<F, coeffType>{
  public:
	using T=Eigen::Matrix<F, Eigen::Dynamic, 1>;
	// Do it like this: We cannot generate EnergyEsPreFetched
	// directly, this has to be done via EnergyEs
	friend EnergyEs<F, coeffType>;

	T nabla(State<coeffType> const &input) const;

	// Allow for polymorphic copy
	EnergyEsPreFetched* clone() const {return new EnergyEsPreFetched(*this);}

	// For sake of completeness, we specify that this requires connections
	int connectionsRequired() const {return numCons;}
private:
    // Make sure this is not manually constructed, but only via
    // EnergyEs. This way, we cannot attribute the wrong CF to a sampler
	explicit EnergyEsPreFetched(Hamiltonian const &H_, int numCons_):
		EnergyCFBaseClass<F, coeffType>(H_), numCons(numCons_){};

	coeffType evaluate(State<coeffType> const &input) const;

	using EnergyCFBaseClass<F, coeffType>::H;
	using EnergyCFBaseClass<F, coeffType>::energy;
	using EnergyCFBaseClass<F, coeffType>::normalizerCoeff;

	int numCons;
};

}

#endif /* SRC_COSTFUNCTIONS_ENERGYESPREFETCHED_HPP_ */
