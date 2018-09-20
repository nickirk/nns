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
#include "../utilities/MatrixTypeDefine.hpp"

namespace networkVMC{

class State;

/** \class CostFunction
 * \brief Pure abstract base class for cost functions
 *
 * Defines the interface used for functions to optimize. The function f to optimize is a real- or complex-valued
 * function over a vector space, and both the function and its derivative have to be supplied.
 * For complex valued functions, the absolute value is minimized.
 */

class CostFunction{
  public:
	/// Type of the derivative
	virtual ~CostFunction(){};
// Two functions have to be present in a cost function: The function itself (calc) and its derivative
// (nabla)

	/**
	 * \brief Interface function for the derivative of the function to optimize
	 *
	 * \param[in] input Input vector as a State object, containing all vector coefficients and corresponding basis vectors
	 * Computes the vector df/dc where c are the vector coefficients of the input state.
	 */
	virtual paraVector nabla(State const &input) const = 0;

	/**
	 * \brief Interface function for the function f to optimize
	 *
	 * \param[in] input Input vector as a State object, containing all vector coefficients and corresponding basis states
	 * Computes f(input)
	 */
	virtual coeffType calc(State const &input) const = 0;

	// Make the CostFunction clonable - as we have multiple inheritance, CRTP is
	// not such a good idea

	/**
	 * \brief Virtual constructor for CostFunction
	 *
	 * \return pointer to a new instance of a CostFunction which is a copy of *this (with the same derived type)
	 */
	virtual CostFunction* clone() const = 0;

	// (formally returns the number of such connections, which is 0 by default)
	/**
	 * \brief Returns the number of additional basis vectors per supplied vector entry required in the input of f
	 *
	 * The input vector typically only has nonzero coefficients with a very few basis vectors. Often, it is helpful to
	 * supply additional, auxiliary entries that can be used to compute the function value. This member function returns
	 * the number of requested auxiliary entries per main entry. Default is 0.
	 */
	virtual int connectionsRequired() const {return 0;}

	// By default, the setup just returns a copy of the CF, so do nothing
	/**
	 * \brief Initializes a CostFunction specialization
	 *
	 * \param[in] sT Type of the Sampler used to create the input vectors
	 * Specializes the CostFunction for the given SamplerType. Exact behavior depends on the used CostFunction
	 * implementation, but this ensures the CostFunction is always suited for a given SamplerType.
	 */
	virtual void setUpCF(SamplerType const &sT) {};
};

}

#endif
