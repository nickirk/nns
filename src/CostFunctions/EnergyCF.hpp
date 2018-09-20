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
 * The function implemented here is the energy expectation value of a supplied Hamiltonian
 */

class EnergyCF: public EnergyCFBaseClass{
  public:
	/**
	 * \param[in] H_ Hamiltonian defining the energy functional
	 */
	explicit EnergyCF(Hamiltonian const &H_):EnergyCFBaseClass(H_){};
	virtual ~EnergyCF(){};
// implementation of the function itself and its derivative
	virtual paraVector nabla(State const &input) const;

	// Allow for polymorphic copy
	virtual EnergyCF* clone() const {return new EnergyCF(*this);}
private:
	/**
	 * \brief Compute the energy expectation value and cache it
	 *
	 * \param[in] input Vector represented by a State object
	 * \return Expectation value of the energy
	 */
	coeffType evaluate(State const &input) const;
};

}

#endif /* SRC_COSTFUNCTIONS_ENERGYCF_HPP_ */
