/*
 * CostFunction.cxx
 *
 *  Created on: Oct 25, 2017
 *      Author: Ke Liao, Kai Guther
 */

#ifndef COST_FUNCTION_DEFINED
#define COST_FUNCTION_DEFINED

#include "State.hpp"
#include <Eigen/Dense>
#include "TypeDefine.hpp"

// Abstract base class for cost functions

class CostFunction{
public:
	CostFunction(){};
	virtual ~CostFunction(){};
// Two functions have to be present in a cost function: The function itself (calc) and its derivative
// (nabla)
	virtual std::vector<Eigen::VectorXd > nabla(std::vector<State> const &input) const = 0;
	virtual double calc(std::vector<State> const &input) const = 0;
};

#endif
