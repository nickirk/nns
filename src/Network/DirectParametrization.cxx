/*
 * DirectParametrization.cxx
 *
 *  Created on: May 2, 2018
 *      Author: guther
 */

#include "DirectParametrization.hpp"
#include "../utilities/State.hpp"

namespace networkVMC {

template <typename T>
DirectParametrization<T>::~DirectParametrization() {
}

template <typename T>
T DirectParametrization<T>::calcNablaPars(State const &inputState,
		nablaType const &dEdC){
  // In the direct parametrization, the inner derivative is the unity matrix
  // so we return just dEdC (the real part)
  int numDets = inputState.size();
  T dEdPars = Eigen::VectorXd::Zero(coeffs.size());
  int j = 0;
  for(int i=0; i<numDets; ++i){
	  // get the index corresponding to the i-th det of inputState
	  j = fullBasis->getIndexByDet(inputState.det(i));
	  // take the first entry of the nabla-type (this is the real part)
	  dEdPars[j] = dEdC[i].real();
  }
  // and return
  return dEdPars;
}

template <typename T>
Eigen::VectorXcd DirectParametrization<T>::getDeriv(detType const &det) const {
  Eigen::VectorXcd dCdW=Eigen::VectorXcd::Zero(coeffs.size()); 
  dCdW(fullBasis->getIndexByDet(det))=coeffType(1.0,0.);
  return dCdW;
}

template <typename T>
Eigen::VectorXcd DirectParametrization<T>::getMarkovDeriv(detType const &det) const {
  Eigen::VectorXcd dCdW=Eigen::VectorXcd::Zero(coeffs.size()); 
  dCdW(fullBasis->getIndexByDet(det))=1./coeffs[fullBasis->getIndexByDet(det)];
  return dCdW;
}
//---------------------------------------------------------------------------//
//instantiate class
//template class DirectParametrization<VecType>;
template class DirectParametrization<VecCType>;
} /* namespace networkVMC */
