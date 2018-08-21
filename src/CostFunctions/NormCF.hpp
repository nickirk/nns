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

////class State;

// Cost function that measures the distance ||psi - psi_0||^2
template <typename F=std::complex<double>, typename coeffType=std::complex<double>>
class NormCF: public CostFunction<F, coeffType>{
  public:
	using T=Eigen::Matrix<F, Eigen::Dynamic, 1>;
	// NormCF construction/destruction
	explicit NormCF(State<coeffType> const &psi_):psi(psi_){};
	virtual ~NormCF(){};
// derivative of ||psi - psi_0||^2 with respect to the coefficients of psi
	T nabla(State<coeffType> const &input) const;
// value of ||psi - psi_0||^2
	coeffType calc(State<coeffType> const &input) const;

	// Allow for polymorphic copy
	virtual CostFunction<F, coeffType>* clone() const {return new NormCF<F, coeffType>(*this);}
private:
	// The reference state psi_0
	State<coeffType> psi;
};

}

#endif /* SRC_COSTFUNCTIONS_NORMCF_HPP_ */
