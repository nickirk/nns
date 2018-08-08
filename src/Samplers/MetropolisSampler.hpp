/*
 * MarkovSampler.hpp
 *
 *  Created on: Jan 23, 2018
 *      Author: guther
 */

#ifndef SRC_SAMPLERS_METROPOLISSAMPLER_HPP_
#define SRC_SAMPLERS_METROPOLISSAMPLER_HPP_

#include "Sampler.hpp"
#include "../Network/Parametrization.hpp"

// Class for sampling determinants using the metropolis algorithm

namespace networkVMC{
template <typename T=VecType>
class MetropolisSampler: public Sampler {
public:
	MetropolisSampler(ExcitationGenerator const &eG_, detType const &HF,
			          Basis const &fullBasis_, Parametrization<T> const &para_, 
                int numDets_ = 100):Sampler(eG_,HF,numDets_),para(&para_),
                lastCoeff(para_.getCoeff(cDet)),fullBasis(&fullBasis_){};
	MetropolisSampler(Hamiltonian const &H_, detType const &HF, 
                Basis const &fullBasis_,Parametrization<T> const &para_, 
                int numDets_ = 100):Sampler(H_,HF,numDets_),para(&para_),
                fullBasis(&fullBasis_),lastCoeff(para_.getCoeff(cDet)){};

  MetropolisSampler(MetropolisSampler const &source):
    Sampler(static_cast<Sampler const&>(source)){
    para = source.para;
    fullBasis = source.fullBasis;
    setReference(getRandomDeterminant(*fullBasis));
  };

  MetropolisSampler(MetropolisSampler &&source):
    Sampler(static_cast<Sampler &>(source)){ 
    swap(*this,source);
    setReference(getRandomDeterminant(*fullBasis));
  };

  friend void swap(MetropolisSampler &a, MetropolisSampler &b){
    swap(static_cast<Sampler&>(a),static_cast<Sampler&>(b));
    std::swap(a.para,b.para);
    std::swap(a.fullBasis, b.fullBasis);
  }

  MetropolisSampler operator=(MetropolisSampler source){
    swap(*this,source);
  // cDet to random
    setReference(getRandomDeterminant(*fullBasis));
    return *this;
  }
	virtual ~MetropolisSampler();
	// create a dynamic polymorphic copy
	virtual MetropolisSampler* clone() const {return new MetropolisSampler(*this);}
	//Do a markov step
	virtual void iterate(coeffType &cI, detType &dI, double &weight, int i);
	// this is a markov-type sampler
	SamplerType type() const {return Markov;}

  void setReference(detType const &a) {
    Sampler::setReference(a);
    lastCoeff = para->getCoeff(a);
  }
	void resetSpecs(){lastCoeff = para->getCoeff(cDet);}
private:
  // sampling depends on the coefficients, as they have to be given alongside the determinants
  Parametrization<T> const *para;
	coeffType lastCoeff;
  Basis const *fullBasis;
};

}
#endif /* SRC_SAMPLERS_METROPOLISSAMPLER_HPP_ */
