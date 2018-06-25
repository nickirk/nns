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

// This is for increased safety, we can set up the CF depending
// on the type of the sampler, to avoid using invalid combinations
	virtual CostFunction const& setUpCF(SamplerType const &sT) const{
		return *this;
	}
};

}

#endif
