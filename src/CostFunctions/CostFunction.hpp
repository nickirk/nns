/*
 * CostFunction.cxx
 *
 *  Created on: Oct 25, 2017
 *      Author: Ke Liao, Kai Guther
 */

#ifndef COST_FUNCTION_DEFINED
#define COST_FUNCTION_DEFINED

#include <Eigen/Dense>
#include "../utilities/TypeDefine.hpp"
#include <memory>

namespace networkVMC{

// Abstract base class for cost functions

class State;

class CostFunction{
public:
	CostFunction(){};
	virtual ~CostFunction(){};
// Two functions have to be present in a cost function: The function itself (calc) and its derivative
// (nabla)
	virtual nablaType nabla(State const &input) const = 0;
	virtual double calc(State const &input) const = 0;

	// Make the CostFunction clonable - as we have multiple inheritance, CRTP is
	// not such a good idea
	virtual CostFunction* clone() const = 0;

	// By default, the setup just returns a copy of the CF
	virtual std::unique_ptr<CostFunction> setUpCF(SamplerType const &sT) const {
		return std::unique_ptr<CostFunction>(this->clone());
	}
};

}

#endif
