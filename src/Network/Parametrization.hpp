/*
 * Parametrization.hpp
 *
 *  Created on: Apr 30, 2018
 *      Author: guther
 */

#ifndef SRC_NETWORK_PARAMETRIZATION_HPP_
#define SRC_NETWORK_PARAMETRIZATION_HPP_

#include "../utilities/TypeDefine.hpp"
#include <Eigen/Dense>
#include "../utilities/StateForward.hpp"
#include "../utilities/Errors.hpp"

// set the default template arguments
#include "ParametrizationForward.hpp"

namespace networkVMC {

// This is the base class for any parametrization of the wave function
// It is an abstract parametrization that can be updated to be optimized
// with respect to a given cost function
template <typename F, typename coeffType>
class Parametrization {
  public:
	Parametrization() = default;
	using T=Eigen::Matrix<F, Eigen::Dynamic ,1>;
    // default the virtual destructor
    virtual ~Parametrization() = default;
    // default the move/copy constructors and assignment operators (not automatically done
    // due to the virtual destructor)
    Parametrization(Parametrization<F, coeffType> &&) = default;
    Parametrization(Parametrization<F, coeffType> const&) = default;
    Parametrization<F, coeffType>& operator=(Parametrization<F, coeffType> &&) = default;
    Parametrization<F, coeffType>& operator=(Parametrization<F, coeffType> const&) = default;

    // It needs to be able to return coefficients somehow
    virtual coeffType getCoeff(detType const &det) const=0; // Can throw an invalidDeterminantError
    // base method for returning the parameters (as a single vector)
    virtual T const& pars() const = 0;
    // We also need write access to the parameters (this is the non-const variant)
    virtual T& pars(){return const_cast<T&>
      (static_cast<Parametrization<F, coeffType> const&>(*this).pars());};

    void writeParsToFile(std::string file) const;
    void readParsFromFile(std::string file);

    int getNumPars(){ return pars().size();};

  // The following function is used when use ListGen or fullSampler
    virtual T calcNablaParsConnected(State<coeffType> const &inputState,
    		T const &dEdC);

  // And this one is for Metropolis or Umbrella sampling
    virtual T calcNablaParsMarkovConnected(State<coeffType> const &inputState,
            T const& dEdC, F const& energy);

    // virtual construction
    virtual Parametrization<F, coeffType>* clone() const = 0;
    virtual Parametrization<F, coeffType>* move_clone() = 0;

    //friend class TrialWfPara<F, coeffType>;
    virtual T getDeriv(detType const &det) const=0;
    //   Obtain the inner derivative dX/dPars with given dX/dC (C are coefficients)
    virtual T getMarkovDeriv(detType const &det) const{
      	  return getDeriv(det)/getCoeff(det);
        };
};

// implement the polymorphic copy for the derived classes
template <typename T, typename coeffType, typename Base>
class ClonableParametrization: public Parametrization<T, coeffType>{
  public:
	// inherit the constructor
	using Parametrization<T, coeffType>::Parametrization;
	virtual ~ClonableParametrization() = default;

	// copy-clone (does not change *this)
	virtual Parametrization<T, coeffType>* clone() const{
		return new Base{static_cast<Base const&>(*this)};
	}

	// move-clone (moves the data to the return pointer)
	virtual Parametrization<T, coeffType>* move_clone(){
		return new Base(std::move(static_cast<Base &>(*this) ) );
	}
};

} /* namespace networkVMC */

#endif /* SRC_NETWORK_PARAMETRIZATION_HPP_ */
