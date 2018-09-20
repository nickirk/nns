/*
 * DirectParametrization.hpp
 *
 *  Created on: May 2, 2018
 *      Author: guther
 */

#ifndef SRC_NETWORK_DIRECTPARAMETRIZATION_HPP_
#define SRC_NETWORK_DIRECTPARAMETRIZATION_HPP_

#include "Parametrization.hpp"
#include "../utilities/TypeDefine.hpp"
#include "../utilities/MatrixTypeDefine.hpp"
#include "../utilities/Errors.hpp"
#include "../HilbertSpace/Basis.hpp"

namespace networkVMC {

// here we parametrize the wave function by its coefficients
// it is mainly for testing purposes
class DirectParametrization: public ClonableParametrization<DirectParametrization> {
  public:
  DirectParametrization(Basis const &fullBasis_):fullBasis(&fullBasis_){
    // start with a random vector
    coeffs = paraVector::Random(fullBasis->size());
  }

  // The coefficient of a determinant is just its entry
  coeffType getCoeff(detType const &det) const{
	  try{
	    auto i = fullBasis->getIndexByDet(det);
	    return coeffs[i];
	  }
	  catch (errors::OutOfRangeError const&){
		  // If the fed determinant is not valid, this is a problem
		  throw errors::InvalidDeterminantError(det);
	  }
  }

  paraVector const& pars() const{
	  return coeffs;
  }
  // inner derivative implementation
  paraVector calcNablaPars(State const &inputState, paraVector const &dEdC);
  paraVector getDeriv(detType const &det) const;
  paraVector getMarkovDeriv(detType const &det) const;

private:
  // The basis to which we refer
  Basis const *fullBasis;
  // directly store the coefficients
  paraVector coeffs;
};

} /* namespace networkVMC */

#endif /* SRC_NETWORK_DIRECTPARAMETRIZATION_HPP_ */
