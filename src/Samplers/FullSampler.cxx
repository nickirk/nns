/*
 * FullSampler.cxx
 *
 *  Created on: May 3, 2018
 *      Author: guther
 */

#include "FullSampler.hpp"

namespace networkVMC {

FullSampler::FullSampler(Hamiltonian const &H_, Basis const &fullBasis_,
		detType const &HF, Parametrization const & para_):
	Sampler(H_,fullBasis_,HF,para_,fullBasis_.getSize()),pos(0){
}

FullSampler::~FullSampler() {
	// TODO Auto-generated destructor stub
}

void FullSampler::iterate(coeffType &cI, detType &dI) const{
  // No sampling involved: We take the next determinant from the list
  dI = fullBasis->getDetByIndex(pos);
  cI = para->getCoeff(dI);
  pos +=1;
  // And in case we reached the end, start anew
  if(pos>=numDets) pos -= numDets;
}

} /* namespace networkVMC */
