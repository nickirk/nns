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
namespace networkVMC{

/**
 * \class NormCF
 * \brief Cost function that measures the L2-distance to a target vector
 *
 * This class implements a CostFunction using the L2-distance to a supplied target vector, optimizing it tries
 * to create a parametrization of the target vector
 */
class NormCF: public CostFunction{
  public:
	// NormCF construction/destruction
	/**
	 * \param psi_ State object representing the target vector
	 */
	explicit NormCF(State const &psi_):psi(psi_){};
	virtual ~NormCF(){};
// derivative of ||psi - psi_0||^2 with respect to the coefficients of psi
	paraVector nabla(State const &input) const;
// value of ||psi - psi_0||^2
	coeffType calc(State const &input) const;

	// Allow for polymorphic copy
	virtual CostFunction* clone() const {return new NormCF(*this);}
private:
	/// The target state
	State psi;
};

}

#endif /* SRC_COSTFUNCTIONS_NORMCF_HPP_ */
