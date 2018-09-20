/*
 * DirectParametrization.cxx
 *
 *  Created on: May 2, 2018
 *      Author: guther
 */

#include "DirectParametrization.hpp"
#include "../utilities/State.hpp"

namespace networkVMC {

paraVector DirectParametrization::getDeriv(detType const &det) const {
  paraVector dCdW=paraVector::Zero(coeffs.size()); 
  dCdW(fullBasis->getIndexByDet(det))=paraType(1.0);
  return dCdW;
}

//---------------------------------------------------------------------------//

paraVector DirectParametrization::getMarkovDeriv(detType const &det) const {
  paraVector dCdW=paraVector::Zero(coeffs.size()); 
  dCdW(fullBasis->getIndexByDet(det))=1./coeffs[fullBasis->getIndexByDet(det)];
  return dCdW;
}

} /* namespace networkVMC */
