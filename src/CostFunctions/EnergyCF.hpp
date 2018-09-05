/*
 * EnergyCF.hpp
 *
 *  Created on: Oct 25, 2017
 *      Author: guther
 */

#ifndef SRC_COSTFUNCTIONS_ENERGYCF_HPP_
#define SRC_COSTFUNCTIONS_ENERGYCF_HPP_

#include <vector>
#include "../utilities/TypeDefine.hpp"
#include <Eigen/Dense>
#include "EnergyCFBaseClass.hpp"
// This cost function tries to minimize the energy expectation value

namespace networkVMC{

class Hamiltonian;

/**
 * \class EnergyCF
 * \brief This CostFunction uses the energy expecation value
 * \tparam F Type of the parameters to optimize
 * \tparam coeffType Type of the vector coefficients of the input vector
 * The function implemented here is the energy expectation value of a supplied Hamiltonian
 */

template <typename F=std::complex<double>, typename coeffType=std::complex<double>>
class EnergyCF: public EnergyCFBaseClass<F, coeffType>{
  public:
	/// Type of the derivative
	using T=Eigen::Matrix<F, Eigen::Dynamic, 1>;
	/**
	 * \param H_ Hamiltonian defining the energy functional
	 */
	explicit EnergyCF(Hamiltonian const &H_):EnergyCFBaseClass<F, coeffType>(H_){};
	virtual ~EnergyCF(){};
// implementation of the function itself and its derivative
	virtual T nabla(State<coeffType> const &input) const;

	// Allow for polymorphic copy
	virtual EnergyCF<F, coeffType>* clone() const {return new EnergyCF<F, coeffType>(*this);}
private:
	/**
	 * \brief Compute the energy expectation value and cache it
	 *
	 * \param input Vector represented by a State object
	 * \return Expectation value of the energy
	 */
	coeffType evaluate(State<coeffType> const &input) const;
	using EnergyCFBaseClass<F, coeffType>::H;
	using EnergyCFBaseClass<F, coeffType>::energy;
	using EnergyCFBaseClass<F, coeffType>::normalizerCoeff;
};

}

#endif /* SRC_COSTFUNCTIONS_ENERGYCF_HPP_ */
