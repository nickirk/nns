/*
 * EnergyCFBaseClass.hpp
 *
 *  Created on: Jun 26, 2018
 *      Author: guther
 */

#ifndef SRC_COSTFUNCTIONS_ENERGYCFBASECLASS_HPP_
#define SRC_COSTFUNCTIONS_ENERGYCFBASECLASS_HPP_

#include "CostFunction.hpp"

namespace networkVMC {

class Hamiltonian;
////class State;

/**
 * \class EnergyCFBaseClass
 * \brief Base class for energy-based CostFunctions
 *
 * CostFunctions based on Energy expectation values are of this type,
 * they have internal cache for the energy and the norm and
 * evaluate them during calls of nabla(), such that calc() comes at no extra cost
 *
*/
template <typename F=std::complex<double>, typename coeffType=std::complex<double>>
class EnergyCFBaseClass: public CostFunction<F, coeffType> {
  public:
	/// Type of the derivative
	using T=Eigen::Matrix<F, Eigen::Dynamic, 1>;
	/**
	 * \param H_ Hamiltonian used for obtaining the energy
	 */
	EnergyCFBaseClass(Hamiltonian const &H_):H(H_), energy(0.0), 
  normalizerCoeff(0.0){};
	virtual ~EnergyCFBaseClass(){};

	// these are the same in all EnergyEs cost functions
	virtual coeffType calc(State<coeffType> const &input) const {return energy;};
	/**
	 * \brief Auxiliary function to get the norm
	 * \return The norm of the last input state
	 */
	virtual double getNormalizer() const{return normalizerCoeff;};

	// The EnergyEs base class is still abstract, as the rules of how to
	// get the energy & nabla are to be specified
	virtual T nabla(State<coeffType> const &input) const = 0;

	virtual EnergyCFBaseClass *clone() const = 0;
protected:
	/// Hamiltonian, to which the energy refers
	Hamiltonian const& H;
	/// cache variable for intermediate results
	mutable coeffType energy;
	/// cache variable for intermediate results
	mutable double normalizerCoeff;
};

} /* namespace networkVMC */

#endif /* SRC_COSTFUNCTIONS_ENERGYCFBASECLASS_HPP_ */
