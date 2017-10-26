/*
 * CostFunction.cxx
 *
 *  Created on: Oct 25, 2017
 *      Author: Ke Liao, Kai Guther
 */

#ifndef COST_FUNCTION_DEFINED
#define COST_FUNCTION_DEFINED

#include "State.hpp"
#include "CoeffType.hpp"

// Abstract base class for cost functions

class CostFunction{
public:
	CostFunction(){};
	virtual ~CostFunction(){};
// To prevent redundant calls to calc(), all calls to the cost function have to be made via
// an evaluator object
	friend class Evaluator;
private:
// Three functions have to be present in a cost function: The function itself (calc), its derivative
// (nabla) and a getter for the function value (getValue)
	virtual std::vector<coeffType > nabla(State const &input) const = 0;
	virtual void calc(State const &input) const = 0;
	virtual double getValue() const = 0;
};

#endif
