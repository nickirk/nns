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
 * \tparam F Type of the parameters to optimize
 * \tparam coeffType Type of the vector coefficients of the input vector
 * This class implements a CostFunction using the L2-distance to a supplied target vector, optimizing it tries
 * to create a parametrization of the target vector
 */
template <typename F=std::complex<double>, typename coeffType=std::complex<double>>
class NormCF: public CostFunction<F, coeffType>{
  public:
	/// Type of the derivative
	using T=Eigen::Matrix<F, Eigen::Dynamic, 1>;
	// NormCF construction/destruction
	/**
	 * \param psi_ State object representing the target vector
	 */
	explicit NormCF(State<coeffType> const &psi_):psi(psi_){};
	virtual ~NormCF(){};
// derivative of ||psi - psi_0||^2 with respect to the coefficients of psi
	T nabla(State<coeffType> const &input) const;
// value of ||psi - psi_0||^2
	coeffType calc(State<coeffType> const &input) const;

	// Allow for polymorphic copy
	virtual CostFunction<F, coeffType>* clone() const {return new NormCF<F, coeffType>(*this);}
private:
	/// The target state
	State<coeffType> psi;
};

}

#endif /* SRC_COSTFUNCTIONS_NORMCF_HPP_ */
