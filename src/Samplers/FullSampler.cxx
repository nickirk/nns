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

template <typename F, typename coeffType>
FullSampler<F, coeffType>::FullSampler(ExcitationGenerator const &eG_, Basis const &fullBasis_,
		Parametrization<F, coeffType> const & para_):
	Sampler<coeffType>(eG_,fullBasis_.size()),pos(0),para(&para_),fullBasis(&fullBasis_){
}

//---------------------------------------------------------------------------//

template <typename F, typename coeffType>
FullSampler<F, coeffType>::FullSampler(Hamiltonian const &H_, Basis const &fullBasis_,
		Parametrization<F, coeffType> const & para_):
	Sampler<coeffType>(H_,detType(),fullBasis_.size()),pos(0),para(&para_),fullBasis(&fullBasis_){
}

//---------------------------------------------------------------------------//


template <typename F, typename coeffType>
FullSampler<F, coeffType>::~FullSampler() {
	// TODO Auto-generated destructor stub
}

//---------------------------------------------------------------------------//

template <typename F, typename coeffType>
void FullSampler<F, coeffType>::iterate(coeffType &cI, detType &dI, double& weight, int i){
  if(i >= numDets) i = 0;
  // No sampling involved: We take the i-th determinant from the list
  dI = fullBasis->getDetByIndex(i);
  weight = 1;
  cI = para->getCoeff(dI);
}
//instantiate class
template class FullSampler<double, double>;
;
template class FullSampler<std::complex<double>, std::complex<double>>;
} /* namespace networkVMC */
