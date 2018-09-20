/*
 * MarkovSampler.hpp
 *
 *  Created on: Jan 23, 2018
 *      Author: guther
 */

#ifndef SRC_SAMPLERS_METROPOLISSAMPLER_HPP_
#define SRC_SAMPLERS_METROPOLISSAMPLER_HPP_

#include "Sampler.hpp"

// Class for sampling determinants using the metropolis algorithm

namespace networkVMC{

class Parametrization;

class MetropolisSampler: public Sampler {
  public:
	MetropolisSampler(ExcitationGenerator const &eG_, detType const &HF,
			          Basis const &fullBasis_, Parametrization const &para_,
                int numDets_ = 100);
	MetropolisSampler(Hamiltonian const &H_, detType const &HF, 
                Basis const &fullBasis_,Parametrization const &para_,
                int numDets_ = 100);

	// when we copy or assign a metropolis sampler, we reset the starting determinant to a random one
  MetropolisSampler(MetropolisSampler const &source):
    Sampler(static_cast<Sampler const&>(source)){
    para = source.para;
    fullBasis = source.fullBasis;
    setReference(randDet());
  };

  MetropolisSampler(MetropolisSampler &&source):
    Sampler(static_cast<Sampler &>(source)){
    swap(*this,source);
    setReference(randDet());
  };

  friend void swap(MetropolisSampler &a, MetropolisSampler &b){
    swap(static_cast<Sampler&>(a),static_cast<Sampler&>(b));
    std::swap(a.para,b.para);
    std::swap(a.fullBasis, b.fullBasis);
  }

  MetropolisSampler operator=(MetropolisSampler source){
    swap(*this,source);
  // cDet to random
    setReference(randDet());
    return *this;
  }
	~MetropolisSampler();
	// create a dynamic polymorphic copy
	MetropolisSampler* clone() const {return new MetropolisSampler(*this);}
	//Do a markov step
	//template <typename coeffType=std::complex<double>>
	void iterate(coeffType &cI, detType &dI, double &weight, int i);
	// this is a markov-type sampler
	SamplerType type() const {return Markov;}

  void setReference(detType const &a);
private:
  // sampling depends on the coefficients, as they have to be given alongside the determinants
  Parametrization const *para;
  coeffType lastCoeff;
  Basis const *fullBasis;
  detType cDet;


  // some random determinant - we have to change this at some point when we remove the basis
  detType randDet() const;
};

}
#endif /* SRC_SAMPLERS_METROPOLISSAMPLER_HPP_ */
