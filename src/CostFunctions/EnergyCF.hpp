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

template <typename F=std::complex<double>, typename coeffType=std::complex<double>>
class EnergyCF: public EnergyCFBaseClass<F, coeffType>{
  public:
	using T=Eigen::Matrix<F, Eigen::Dynamic, 1>;
	// Here, we need to supply a Hamiltonian
	explicit EnergyCF(Hamiltonian const &H_, int numCons_=20):EnergyCFBaseClass<F, coeffType>(H_),
   numCons(numCons_){};
	virtual ~EnergyCF(){};
// implementation of the function itself and its derivative
	virtual T nabla(State<coeffType> const &input) const;

	// Allow for polymorphic copy
	virtual EnergyCF<F, coeffType>* clone() const {return new EnergyCF<F, coeffType>(*this);}
private:
	coeffType evaluate(State<coeffType> const &input) const;
	using EnergyCFBaseClass<F, coeffType>::H;
	using EnergyCFBaseClass<F, coeffType>::energy;
	using EnergyCFBaseClass<F, coeffType>::normalizerCoeff;
  int numCons;
};

}

#endif /* SRC_COSTFUNCTIONS_ENERGYCF_HPP_ */
