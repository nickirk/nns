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
#include "CostFunction.hpp"
#include "Eigen/Dense"
#include "EnergyCFBaseClass.hpp"
// This cost function tries to minimize the energy expectation value

namespace networkVMC{

class Hamiltonian;

class EnergyCF: public EnergyCFBaseClass{
public:
	// Here, we need to supply a Hamiltonian
	explicit EnergyCF(Hamiltonian const &H_):EnergyCFBaseClass(H_){};
	virtual ~EnergyCF(){};
// implementation of the function itself and its derivative
	virtual nablaType nabla(State const &input) const;

	// Allow for polymorphic copy
	virtual CostFunction* clone() const {return new EnergyCF(*this);}
private:
	double evaluate(State const &input) const;
};

}

#endif /* SRC_COSTFUNCTIONS_ENERGYCF_HPP_ */
