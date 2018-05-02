/*
 * Parametrization.hpp
 *
 *  Created on: Apr 30, 2018
 *      Author: guther
 */

#ifndef SRC_NETWORK_PARAMETRIZATION_HPP_
#define SRC_NETWORK_PARAMETRIZATION_HPP_

#include "../utilities/TypeDefine.hpp"
#include "../utilities/State.hpp"
#include "../CostFunctions/CostFunction.hpp"

namespace networkVMC {

// This is the base class for any parametrization of the wave function
// It is an abstract parametrization that can be updated to be optimized
// with respect to a given cost function

class Parametrization {
public:
  Parametrization(){};
  virtual ~Parametrization(){};
  // We need access to the parameters (this is the non-const variant)
  virtual VecType& pars(){return const_cast<VecType>
    (static_cast<Parametrization const&>(*this).pars());};
 // It needs to be able to return coefficients somehow
  virtual coeffType getCoeff(detType const &det) const=0;
  // base method for returning the parameters (as a single vector)
  virtual VecType const& pars() const = 0;
// Obtain the inner derivative dX/dPars with given dX/dC (C are coefficients)
  virtual VecType calcNablaPars(
		  State const &input,
		  nablaType const &outerDerivative) = 0;
};

} /* namespace networkVMC */

#endif /* SRC_NETWORK_PARAMETRIZATION_HPP_ */
