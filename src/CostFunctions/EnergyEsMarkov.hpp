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
#include "EnergyCFBaseClass.hpp"

#include "EnergyEsForward.hpp"

namespace networkVMC{

class Hamiltonian;

/**
 * \class EnergyEsMarkov
 *
 * \brief Implementation of the EnergyEs for Markov-Chain sampling
 *
 * \tparam F Type of the parameters to optimize
 * \tparam coeffType Type of the vector coefficients of the input vector
 *
 * This class implements the EnergyEs CostFunction (expectation value with respect to stochastic sampling)
 * for a Markov-Chain sampling
 */
template <typename F=std::complex<double>, typename coeffType=std::complex<double>>
class EnergyEsMarkov: public EnergyCFBaseClass<F, coeffType>{
  public:
	/// Type of the derivative
	using T=Eigen::Matrix<F, Eigen::Dynamic, 1>;
	/// EnergyEsMarkov can only be generated by EnergyEs
	friend EnergyEs<F, coeffType>;

	coeffType calc(State<coeffType> const &input) const {return energy;};
	T nabla(State<coeffType> const &input) const;

	// Allow for polymorphic copy
	EnergyEsMarkov* clone() const {return new EnergyEsMarkov(*this);}

	// For sake of completeness, we specify that this requires connections
	int connectionsRequired() const {return numCons;}
private:
    // Make sure this is not manually constructed, but only via
    // EnergyEs. This way, we cannot attribute the wrong CF to a sampler
	/**
	 * \param[in] H_ Hamiltonian used to obtain the energy
	 * \param[in] numCons_ Number of matrix elements taken into account per basis vector
	 * Can only be called by the friend class
	 */
	explicit EnergyEsMarkov(Hamiltonian const &H_,int numCons_):
		EnergyCFBaseClass<F, coeffType>(H_),numCons(numCons_){};

	/**
	 * \brief Compute the energy expectation value and cache it
	 *
	 * \param[in] input Vector represented by a State object
	 * \return Expectation value of the energy
	 */
	coeffType evaluate(State<coeffType> const &input) const;

	using EnergyCFBaseClass<F, coeffType>::H;
	using EnergyCFBaseClass<F, coeffType>::energy;
	using EnergyCFBaseClass<F, coeffType>::normalizerCoeff;
	/// Number of matrix elements to sample per basis vector
	int numCons;
};

}

#endif /* SRC_COSTFUNCTIONS_ENERGYESMARKOV_HPP_ */
