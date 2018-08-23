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
#include "../Network/ParametrizationForward.hpp"
#include "../Network/TrialWfParaForward.hpp"

namespace networkVMC{

template <typename F=std::complex<double>, typename coeffType=std::complex<double>>
class UmbrellaSampler: public Sampler<coeffType> {
  public:
	UmbrellaSampler(ExcitationGenerator const &eG_, detType const &HF,
			    Basis const &fullBasis_, TrialWfPara<F, coeffType> const &para_,
                int numDets_ = 100);
	UmbrellaSampler(Hamiltonian const &H_, detType const &HF,
                Basis const &fullBasis_,TrialWfPara<F, coeffType> const &para_,
                int numDets_ = 100);

    void iterate(coeffType &cI, detType &dI, double &weight, int i);

    Sampler<coeffType>* clone() const {return new UmbrellaSampler<F, coeffType>(*this);}

    void setReference(detType const &start){internalSampler.setReference(start);}

    SamplerType type() const{return internalSampler.type();}
  private:
    using Sampler<coeffType>::numDets;
    detType cDet;
    MetropolisSampler<F, coeffType> internalSampler;
    TrialWfPara<F,coeffType> const *para;
};

} //for namespace

#endif /* SRC_SAMPLERS_UMBRELLASAMPLER_HPP_ */
