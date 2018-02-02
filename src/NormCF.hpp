/*
 * NormCF.hpp
 *
 *  Created on: Oct 26, 2017
 *      Author: guther
 */

#ifndef SRC_NORMCF_HPP_
#define SRC_NORMCF_HPP_

#include "CostFunction.hpp"
#include "State.hpp"
#include "CoeffType.hpp"

// Cost function that measures the distance ||psi - psi_0||^2

class NormCF: public CostFunction{
public:
	explicit NormCF(std::vector<State> const &psi_):psi(psi_){};
	std::vector<Eigen::VectorXd > nabla(std::vector<State> const &input) const;
	double calc(std::vector<State> const &input) const;
	double getValue(std::vector<State> const &input) const{return calc(input);}
private:
	std::vector<State> psi;
};


#endif /* SRC_NORMCF_HPP_ */
