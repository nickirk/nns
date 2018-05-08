/*
 * Parametrization.hpp
 *
 *  Created on: Apr 30, 2018
 *      Author: guther
 */

#ifndef SRC_NETWORK_PARAMETRIZATION_HPP_
#define SRC_NETWORK_PARAMETRIZATION_HPP_

#include "../utilities/TypeDefine.hpp"

namespace networkVMC {

// Forward declaration for more efficient compilation
class State;

// This is the base class for any parametrization of the wave function
// It is an abstract parametrization that can be updated to be optimized
// with respect to a given cost function

class Parametrization {
public:
  Parametrization(){};
  virtual ~Parametrization(){};
 // It needs to be able to return coefficients somehow
  virtual coeffType getCoeff(detType const &det) const=0; // Can throw an invalidDeterminantError
  // base method for returning the parameters (as a single vector)
  virtual VecType const& pars() const = 0;
  // We also need write access to the parameters (this is the non-const variant)
  virtual VecType& pars(){return const_cast<VecType&>
    (static_cast<Parametrization const&>(*this).pars());};
// Obtain the inner derivative dX/dPars with given dX/dC (C are coefficients)
  virtual VecType calcNablaPars(
		  State const &input,
		  nablaType const &outerDerivative) = 0;

// The following features are experimental and not essential to the interface
  // Some other derivative
  virtual VecType calcNablaParsConnected(State const &inputState, nablaType const& dEdC)
  	  {return calcNablaPars(inputState,dEdC);}

  // stochastic reconfiguration derivative
  virtual Eigen::MatrixXcd calcdCdwSR(
    State const &outputState
  ){};

};

} /* namespace networkVMC */

#endif /* SRC_NETWORK_PARAMETRIZATION_HPP_ */
