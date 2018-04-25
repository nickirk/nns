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

// Cost function that measures the distance ||psi - psi_0||^2

class NormCF: public CostFunction{
public:
	explicit NormCF(std::vector<State> const &psi_):psi(psi_){};
// derivative of ||psi - psi_0||^2 with respect to the coefficients of psi
	std::vector<Eigen::VectorXd > nabla(std::vector<State> const &input) const;
// value of ||psi - psi_0||^2
	double calc(std::vector<State> const &input) const;
	double getValue(std::vector<State> const &input) const{return calc(input);}
private:
	std::vector<State> psi;
};


#endif /* SRC_COSTFUNCTIONS_NORMCF_HPP_ */
