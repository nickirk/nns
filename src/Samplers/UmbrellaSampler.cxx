/*
 * UmbrellaSampler.cxx
 *
 *  Created on: Aug 22, 2018
 *      Author: liao
 */
#include "UmbrellaSampler.hpp"
#include "../Network/TrialWfPara.hpp"

namespace networkVMC{

template <typename F, typename coeffType>
UmbrellaSampler<F, coeffType>::UmbrellaSampler(ExcitationGenerator const &eG_, detType const &HF,
			    Basis const &fullBasis_, TrialWfPara<F, coeffType> const &para_,
                int numDets_):Sampler<coeffType>(eG_,numDets_),cDet(HF),
				internalSampler(eG_, HF, fullBasis_, para_.base(), numDets_),para(&para_){};

template <typename F, typename coeffType>
UmbrellaSampler<F, coeffType>::UmbrellaSampler(Hamiltonian const &H_, detType const &HF,
                Basis const &fullBasis_,TrialWfPara<F, coeffType> const &para_,
                int numDets_ ):Sampler<coeffType>(H_, HF, numDets_),cDet(HF),
				internalSampler(H_,HF, fullBasis_, para_.base(), numDets_),
				para(&para_){};

//---------------------------------------------------------------------------//

template<typename F, typename coeffType>
void UmbrellaSampler<F,coeffType>::iterate(coeffType &cI, detType &dI, double &weight,
	    int i){
	internalSampler.iterate(cI,dI,weight,i);
	auto tI= para->getTrialCoeff(dI);
	cI *= tI;
	weight *= std::norm(tI);
}

// instantiate it
template class UmbrellaSampler<double, double>;
template class UmbrellaSampler<std::complex<double>, std::complex<double>>;

}
