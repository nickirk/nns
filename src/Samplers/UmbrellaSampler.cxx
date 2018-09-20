/*
 * UmbrellaSampler.cxx
 *
 *  Created on: Aug 22, 2018
 *      Author: liao
 */
#include "UmbrellaSampler.hpp"
#include "../Network/TrialWfPara.hpp"

namespace networkVMC{

UmbrellaSampler::UmbrellaSampler(ExcitationGenerator const &eG_, detType const &HF,
			    Basis const &fullBasis_, TrialWfPara const &para_,
                int numDets_):Sampler(eG_,numDets_),cDet(HF),
				internalSampler(eG_, HF, fullBasis_, para_.base(), numDets_),para(&para_){};

UmbrellaSampler::UmbrellaSampler(Hamiltonian const &H_, detType const &HF,
                Basis const &fullBasis_,TrialWfPara const &para_,
                int numDets_ ):Sampler(H_, HF, numDets_),cDet(HF),
				internalSampler(H_,HF, fullBasis_, para_.base(), numDets_),
				para(&para_){};

//---------------------------------------------------------------------------//

void UmbrellaSampler::iterate(coeffType &cI, detType &dI, double &weight,
	    int i){
	internalSampler.iterate(cI,dI,weight,i);
	auto tI= para->getTrialCoeff(dI);
	cI *= tI;
	weight *= std::norm(tI);
}


}
