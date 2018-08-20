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
#include "../utilities/Errors.hpp"
#include "../HilbertSpace/Basis.hpp"

namespace networkVMC {


// here we parametrize the wave function by its coefficients
// it is mainly for testing purposes
template <typename T=VecType>
class DirectParametrization: public ClonableParametrization<T,DirectParametrization<T> > {
public:
  DirectParametrization(Basis const &fullBasis_):fullBasis(&fullBasis_){
    // start with a random vector
    coeffs = Eigen::VectorXd::Random(fullBasis->size());
  }

  // The coefficient of a determinant is just its entry
  coeffType getCoeff(detType const &det) const{
	  try{
	    auto i = fullBasis->getIndexByDet(det);
	    return coeffs[i];
	  }
	  catch (OutOfRangeError const&){
		  // If the fed determinant is not valid, this is a problem
		  throw InvalidDeterminantError(det);
	  }
  }

  T const& pars() const{
	  return coeffs;
  }
  // inner derivative implementation
  T calcNablaPars(State const &inputState, nablaType const &dEdC);
  Eigen::VectorXcd getDeriv(detType const &det) const;
  Eigen::VectorXcd getMarkovDeriv(detType const &det) const;
  // Some other derivative
  // T calcNablaParsConnected(State const &inputState, nablaType const& dEdC){};
  // stochastic reconfiguration derivative
  // Eigen::MatrixXcd calcdCdwSR(State const &outputState){};

private:
  // The basis to which we refer
  Basis const *fullBasis;
  // directly store the coefficients
  T coeffs;
};

} /* namespace networkVMC */

#endif /* SRC_NETWORK_DIRECTPARAMETRIZATION_HPP_ */
