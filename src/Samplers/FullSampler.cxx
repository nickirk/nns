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
template <typename T>
FullSampler<T>::FullSampler(ExcitationGenerator const &eG_, Basis const &fullBasis_,
		Parametrization<T> const & para_):
	Sampler(eG_,fullBasis_,detType(),fullBasis_.getSize()),pos(0),para(&para_){
}

//---------------------------------------------------------------------------//

template <typename T>
FullSampler<T>::FullSampler(Hamiltonian const &H_, Basis const &fullBasis_,
		Parametrization<T> const & para_):
	Sampler(H_,fullBasis_,detType(),fullBasis_.getSize()),pos(0),para(&para_){
}

//---------------------------------------------------------------------------//


template <typename T>
FullSampler<T>::~FullSampler() {
	// TODO Auto-generated destructor stub
}

//---------------------------------------------------------------------------//

template <typename T>
void FullSampler<T>::iterate(coeffType &cI, detType &dI, double& weight, int i){
  if(i >= numDets) i = 0;
  // No sampling involved: We take the i-th determinant from the list
  dI = fullBasis->getDetByIndex(i);
  weight = 1;
  cI = para->getCoeff(dI);
}
//instantiate class
template class FullSampler<VecType>;
template class FullSampler<VecCType>;
} /* namespace networkVMC */
