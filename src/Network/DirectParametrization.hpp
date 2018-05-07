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
#include "../HilbertSpace/Basis.hpp"

namespace networkVMC {


// here we parametrize the wave function by its coefficients
// it is mainly for testing purposes

class DirectParametrization: public Parametrization {
public:
  DirectParametrization(Basis const &fullBasis_):fullBasis(fullBasis_){
    // start with a random vector
    coeffs = Eigen::VectorXd::Random(fullBasis.getSize());
  }
  virtual ~DirectParametrization();

  // The coefficient of a determinant is just its entry
  coeffType getCoeff(detType const &det) const{
	  auto i = fullBasis.getIndexByDet(det);
	  return coeffs[i];
  }

  VecType const& pars() const{
	  return coeffs;
  }
  // inner derivative implementation
  VecType calcNablaPars(State const &inputState, nablaType const &dEdC);
  // Some other derivative
   VecType calcNablaParsConnected(State const &inputState, nablaType const& dEdC){};
   // stochastic reconfiguration derivative
   Eigen::MatrixXcd calcdCdwSR(State const &outputState){};

private:
  Basis const &fullBasis;
  VecType coeffs;
};

} /* namespace networkVMC */

#endif /* SRC_NETWORK_DIRECTPARAMETRIZATION_HPP_ */
