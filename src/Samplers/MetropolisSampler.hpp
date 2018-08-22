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
template <typename F=std::complex<double>, typename coeffType=std::complex<double>>
class MetropolisSampler: public Sampler<coeffType> {
  public:
	MetropolisSampler(ExcitationGenerator const &eG_, detType const &HF,
			          Basis const &fullBasis_, Parametrization<F, coeffType> const &para_,
                int numDets_ = 100);//:Sampler<F, coeffType>(eG_,HF,numDets_),para(&para_),
                //lastCoeff(para_.getCoeff(cDet)),fullBasis(&fullBasis_){};
	MetropolisSampler(Hamiltonian const &H_, detType const &HF, 
                Basis const &fullBasis_,Parametrization<F, coeffType> const &para_,
                int numDets_ = 100);//:Sampler<F, coeffType>(H_,HF,numDets_),para(&para_),
                //fullBasis(&fullBasis_),lastCoeff(para_.getCoeff(cDet)){};

  MetropolisSampler(MetropolisSampler<F, coeffType> const &source):
    Sampler<coeffType>(static_cast<Sampler<coeffType> const&>(source)){
    para = source.para;
    fullBasis = source.fullBasis;
    //setReference(getRandomDeterminant(*fullBasis));
    setReference(source.getDet());
  };

  MetropolisSampler(MetropolisSampler<F, coeffType> &&source):
    Sampler<coeffType>(static_cast<Sampler<coeffType> &>(source)){
    swap(*this,source);
    //setReference(getRandomDeterminant(*fullBasis));
    setReference(source.getDet());
  };

  friend void swap(MetropolisSampler<F, coeffType> &a, MetropolisSampler<F, coeffType> &b){
    swap(static_cast<Sampler<coeffType>&>(a),static_cast<Sampler<coeffType>&>(b));
    std::swap(a.para,b.para);
    std::swap(a.fullBasis, b.fullBasis);
  }

  MetropolisSampler operator=(MetropolisSampler<F, coeffType> source){
    swap(*this,source);
  // cDet to random
    //setReference(getRandomDeterminant(*fullBasis));
    setReference(source.getDet());
    return *this;
  }
	virtual ~MetropolisSampler();
	// create a dynamic polymorphic copy
	virtual MetropolisSampler<F, coeffType>* clone() const {return new MetropolisSampler<F, coeffType>(*this);}
	//Do a markov step
	//template <typename coeffType=std::complex<double>>
	virtual void iterate(coeffType &cI, detType &dI, double &weight, int i);
	// this is a markov-type sampler
	SamplerType type() const {return Markov;}

  void setReference(detType const &a) {
    Sampler<coeffType>::setReference(a);
    lastCoeff = para->getCoeff(a);
  }
private:
  // sampling depends on the coefficients, as they have to be given alongside the determinants
  Parametrization<F, coeffType> const *para;
  coeffType lastCoeff;
  Basis const *fullBasis;
  using Sampler<coeffType>::numDets;
  using Sampler<coeffType>::cDet;
};

}
#endif /* SRC_SAMPLERS_METROPOLISSAMPLER_HPP_ */
