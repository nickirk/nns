/*
 * MarkovSampler.hpp
 *
 *  Created on: Jan 23, 2018
 *      Author: guther
 */

#ifndef SRC_SAMPLERS_METROPOLISSAMPLER_HPP_
#define SRC_SAMPLERS_METROPOLISSAMPLER_HPP_

#include "Sampler.hpp"
#include "../Network/ParametrizationForward.hpp"

// Class for sampling determinants using the metropolis algorithm

namespace networkVMC{
template <typename F=std::complex<double>, typename coeffType=std::complex<double>>
class MetropolisSampler: public Sampler<coeffType> {
  public:
	MetropolisSampler(ExcitationGenerator const &eG_, detType const &HF,
			          Basis const &fullBasis_, Parametrization<F, coeffType> const &para_,
                int numDets_ = 100);
	MetropolisSampler(Hamiltonian const &H_, detType const &HF, 
                Basis const &fullBasis_,Parametrization<F, coeffType> const &para_,
                int numDets_ = 100);

	// when we copy or assign a metropolis sampler, we reset the starting determinant to a random one
  MetropolisSampler(MetropolisSampler<F, coeffType> const &source):
    Sampler<coeffType>(static_cast<Sampler<coeffType> const&>(source)){
    para = source.para;
    fullBasis = source.fullBasis;
    setReference(randDet());
  };

  MetropolisSampler(MetropolisSampler<F, coeffType> &&source):
    Sampler<coeffType>(static_cast<Sampler<coeffType> &>(source)){
    swap(*this,source);
    setReference(randDet());
  };

  friend void swap(MetropolisSampler<F, coeffType> &a, MetropolisSampler<F, coeffType> &b){
    swap(static_cast<Sampler<coeffType>&>(a),static_cast<Sampler<coeffType>&>(b));
    std::swap(a.para,b.para);
    std::swap(a.fullBasis, b.fullBasis);
  }

  MetropolisSampler operator=(MetropolisSampler<F, coeffType> source){
    swap(*this,source);
  // cDet to random
    setReference(randDet());
    return *this;
  }
	~MetropolisSampler();
	// create a dynamic polymorphic copy
	MetropolisSampler<F, coeffType>* clone() const {return new MetropolisSampler<F, coeffType>(*this);}
	//Do a markov step
	//template <typename coeffType=std::complex<double>>
	void iterate(coeffType &cI, detType &dI, double &weight, int i);
	// this is a markov-type sampler
	SamplerType type() const {return Markov;}

  void setReference(detType const &a);
private:
  // sampling depends on the coefficients, as they have to be given alongside the determinants
  Parametrization<F, coeffType> const *para;
  coeffType lastCoeff;
  Basis const *fullBasis;
  using Sampler<coeffType>::numDets;
  detType cDet;


  // some random determinant - we have to change this at some point when we remove the basis
  detType randDet() const;
};

}
#endif /* SRC_SAMPLERS_METROPOLISSAMPLER_HPP_ */
