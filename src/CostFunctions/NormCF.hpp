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
#include "CostFunction.hpp"
#include "../utilities/State.hpp"

namespace networkVMC{

// Cost function that measures the distance ||psi - psi_0||^2

class NormCF: public CostFunction{
public:
	explicit NormCF(State const &psi_):psi(psi_){};
// derivative of ||psi - psi_0||^2 with respect to the coefficients of psi
	nablaType nabla(State const &input) const;
// value of ||psi - psi_0||^2
	double calc(State const &input) const;
private:
	State psi;
};

}

#endif /* SRC_COSTFUNCTIONS_NORMCF_HPP_ */
