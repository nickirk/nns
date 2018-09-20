/*
 * UmbrellaSampler.hpp
 *
 *  Created on: Aug 22, 2018
 *      Author: liao
 */

#ifndef SRC_SAMPLERS_UMBRELLASAMPLER_HPP_
#define SRC_SAMPLERS_UMBRELLASAMPLER_HPP_

#include "Sampler.hpp"
#include "MetropolisSampler.hpp"

namespace networkVMC{

class TrialWfPara;

class UmbrellaSampler: public Sampler {
  public:
	UmbrellaSampler(ExcitationGenerator const &eG_, detType const &HF,
			    Basis const &fullBasis_, TrialWfPara const &para_,
                int numDets_ = 100);
	UmbrellaSampler(Hamiltonian const &H_, detType const &HF,
                Basis const &fullBasis_,TrialWfPara const &para_,
                int numDets_ = 100);

    void iterate(coeffType &cI, detType &dI, double &weight, int i);

    Sampler* clone() const {return new UmbrellaSampler(*this);}

    void setReference(detType const &start){internalSampler.setReference(start);}

    SamplerType type() const{return internalSampler.type();}
  private:
    using Sampler::numDets;
    detType cDet;
    MetropolisSampler internalSampler;
    TrialWfPara const *para;
};

} //for namespace

#endif /* SRC_SAMPLERS_UMBRELLASAMPLER_HPP_ */
