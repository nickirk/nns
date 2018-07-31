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
template<typename T=VecType>
class Parametrization {
public:
  Parametrization(){};
  virtual ~Parametrization(){};
 // It needs to be able to return coefficients somehow
  virtual coeffType getCoeff(detType const &det) const=0; // Can throw an invalidDeterminantError
  // base method for returning the parameters (as a single vector)
  virtual T const& pars() const = 0;
  // We also need write access to the parameters (this is the non-const variant)
  virtual T& pars(){return const_cast<T&>
    (static_cast<Parametrization const&>(*this).pars());};
// Obtain the inner derivative dX/dPars with given dX/dC (C are coefficients)
  virtual T calcNablaPars(
		  State const &input,
		  nablaType const &outerDerivative) = 0;

  virtual Eigen::VectorXcd getMarkovDeriv(detType const &det) const{
	  return Eigen::VectorXcd();
  };
  virtual Eigen::VectorXcd getDeriv(detType const &det) const{
	  return Eigen::VectorXcd();
  };
// The following features are experimental and not essential to the interface
  // Some other derivative
  virtual T calcNablaParsConnected(State const &inputState, nablaType const& dEdC)
  	  {return calcNablaPars(inputState,dEdC);}

  virtual T calcNablaParsMarkovConnected(State const &inputState, nablaType const& dEdC, double const& energy)
  	  {return calcNablaPars(inputState,dEdC);}

  // stochastic reconfiguration derivative
  virtual Eigen::MatrixXcd calcdCdwSR(
    State const &outputState
  ){return Eigen::MatrixXd();};

};


} /* namespace networkVMC */

#endif /* SRC_NETWORK_PARAMETRIZATION_HPP_ */
