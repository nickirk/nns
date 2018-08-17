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
  // default the big five - required for virtual destructor
  virtual ~Parametrization() = default;

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
		  nablaType const &outerDerivative) = 0; // can throw a SizeMismatchError if input and outerDerivative
  	  	  	  	  	  	  	  	  	  	  	  	 // have different size

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

  // virtual construction
  virtual Parametrization<T>* clone() const = 0;
  virtual Parametrization<T>* move_clone() = 0;

};

// implement the polymorphic copy for the derived classes
template<typename T, typename Base>
class ClonableParametrization: public Parametrization<T>{
public:
	// inherit the constructor
	using Parametrization<T>::Parametrization;
	virtual ~ClonableParametrization() = default;

	// copy-clone (does not change *this)
	virtual Parametrization<T>* clone() const{
		return new Base{static_cast<Base const&>(*this)};
	}

	// move-clone (moves the data to the return pointer)
	virtual Parametrization<T>* move_clone(){
		return new Base(std::move(static_cast<Base &>(*this) ) );
	}
};


} /* namespace networkVMC */

#endif /* SRC_NETWORK_PARAMETRIZATION_HPP_ */
