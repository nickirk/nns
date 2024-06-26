/*
 * FullSampler.cxx
 *
 *  Created on: May 3, 2018
 *      Author: guther
 */

#include "FullSampler.hpp"
#include "../HilbertSpace/Basis.hpp"
#include "../Network/Parametrization.hpp"

namespace networkVMC {

FullSampler::FullSampler(ExcitationGenerator const &eG_, Basis const &fullBasis_,
		Parametrization const & para_):
	Sampler(eG_,fullBasis_.size()),pos(0),para(&para_),fullBasis(&fullBasis_){
}

//---------------------------------------------------------------------------//

FullSampler::FullSampler(Hamiltonian const &H_, Basis const &fullBasis_,
		Parametrization const & para_):
	Sampler(H_,detType(),fullBasis_.size()),pos(0),para(&para_),fullBasis(&fullBasis_){
}

//---------------------------------------------------------------------------//


FullSampler::~FullSampler() {
	// TODO Auto-generated destructor stub
}

//---------------------------------------------------------------------------//

void FullSampler::iterate(coeffType &cI, detType &dI, double& weight, int i){
  if(i >= numDets) i = 0;
  // No sampling involved: We take the i-th determinant from the list
  dI = fullBasis->getDetByIndex(i);
  weight = 1;
  cI = para->getCoeff(dI);
}
} /* namespace networkVMC */
