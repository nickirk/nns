/*
 * DirectParametrization.cxx
 *
 *  Created on: May 2, 2018
 *      Author: guther
 */

#include "DirectParametrization.hpp"
#include "../utilities/State.hpp"

namespace networkVMC {

template <typename F, typename coeffType>
DirectParametrization<F, coeffType>::T DirectParametrization<F, coeffType>::getDeriv(detType const &det) const {
  T dCdW=T::Zero(coeffs.size()); 
  dCdW(fullBasis->getIndexByDet(det))=F(1.0);
  return dCdW;
}

template <typename F, typename coeffType>
DirectParametrization<F, coeffType>::T DirectParametrization<F, coeffType>::getMarkovDeriv(detType const &det) const {
  T dCdW=T::Zero(coeffs.size()); 
  dCdW(fullBasis->getIndexByDet(det))=1./coeffs[fullBasis->getIndexByDet(det)];
  return dCdW;
}
//---------------------------------------------------------------------------//
//instantiate class
template class DirectParametrization<double, double>;
template class DirectParametrization<std::complex<double>, std::complex<double>>;
} /* namespace networkVMC */
