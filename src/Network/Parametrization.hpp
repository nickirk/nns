/*
 * Parametrization.hpp
 *
 *  Created on: Apr 30, 2018
 *      Author: guther
 */

#ifndef SRC_NETWORK_PARAMETRIZATION_HPP_
#define SRC_NETWORK_PARAMETRIZATION_HPP_

#include "../utilities/TypeDefine.hpp"
#include "../utilities/MatrixTypeDefine.hpp"
#include <Eigen/Dense>
#include "../utilities/Errors.hpp"

namespace networkVMC {

class State;

// This is the base class for any parametrization of the wave function
// It is an abstract parametrization that can be updated to be optimized
// with respect to a given cost function
class Parametrization {
  public:
	Parametrization() = default;
    // default the virtual destructor
    virtual ~Parametrization() = default;
    // default the move/copy constructors and assignment operators (not automatically done
    // due to the virtual destructor)
    Parametrization(Parametrization &&) = default;
    Parametrization(Parametrization const&) = default;
    Parametrization& operator=(Parametrization &&) = default;
    Parametrization& operator=(Parametrization const&) = default;

    // It needs to be able to return coefficients somehow
    virtual coeffType getCoeff(detType const &det) const=0; // Can throw an invalidDeterminantError
    // base method for returning the parameters (as a single vector)
    virtual paraVector const& pars() const = 0;
    // We also need write access to the parameters (this is the non-const variant)
    virtual paraVector& pars(){return const_cast<paraVector&>
      (static_cast<Parametrization const&>(*this).pars());};

    void writeParsToFile(std::string file) const;
    void readParsFromFile(std::string file);

    int getNumPars(){ return pars().size();};

  // The following function is used when use ListGen or fullSampler
    virtual paraVector calcNablaParsConnected(State const &inputState,
    		paraVector const &dEdC);

  // And this one is for Metropolis or Umbrella sampling
    virtual paraVector calcNablaParsMarkovConnected(State const &inputState,
            paraVector const& dEdC, paraType const& energy);

    // virtual construction
    virtual Parametrization* clone() const = 0;
    virtual Parametrization* moveClone() = 0;

    //friend class TrialWfPara<F, coeffType>;
    virtual paraVector getDeriv(detType const &det) const=0;
    //   Obtain the inner derivative dX/dPars with given dX/dC (C are coefficients)
    virtual paraVector getMarkovDeriv(detType const &det) const{
      	  return getDeriv(det)/getCoeff(det);
        };
};

// implement the polymorphic copy for the derived classes
template <typename Base>
class ClonableParametrization: public Parametrization{
  public:
	// inherit the constructor
	using Parametrization::Parametrization;
	virtual ~ClonableParametrization() = default;

	// copy-clone (does not change *this)
	virtual Parametrization* clone() const{
		return new Base{static_cast<Base const&>(*this)};
	}

	// move-clone (moves the data to the return pointer)
	virtual Parametrization* moveClone(){
		return new Base(std::move(static_cast<Base &>(*this) ) );
	}
};

} /* namespace networkVMC */

#endif /* SRC_NETWORK_PARAMETRIZATION_HPP_ */
