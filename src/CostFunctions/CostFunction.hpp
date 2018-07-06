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
class CFPars;

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

	// By default, the setup just returns a copy of the CF, so do nothing
	virtual void setUpCF(SamplerType const &sT) {};

	// If we need connected determinants separately contained in the state
	// (formally returns the number of such connections, which is 0 by default)
	virtual int connectionsRequired() const {return 0;}
};

}

#endif
