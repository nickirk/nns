/*
 * NormCF.hpp
 *
 *  Created on: Oct 26, 2017
 *      Author: guther
 */

#ifndef SRC_COSTFUNCTIONS_NORMCF_HPP_
#define SRC_COSTFUNCTIONS_NORMCF_HPP_

#include <Eigen/Dense>
#include "../utilities/TypeDefine.hpp"
#include "../utilities/State.hpp"
#include "CostFunction.hpp"
#include "CostFunction.hpp"

namespace networkVMC{

class State;

// Cost function that measures the distance ||psi - psi_0||^2

class NormCF: public CostFunction{
public:
	// NormCF construction/destruction
	explicit NormCF(State const &psi_):psi(psi_){};
	virtual ~NormCF(){};
// derivative of ||psi - psi_0||^2 with respect to the coefficients of psi
	nablaType nabla(State const &input) const;
// value of ||psi - psi_0||^2
	coeffType calc(State const &input) const;

	// Allow for polymorphic copy
	virtual CostFunction* clone() const {return new NormCF(*this);}
private:
	// The reference state psi_0
	State psi;
};

}

#endif /* SRC_COSTFUNCTIONS_NORMCF_HPP_ */
