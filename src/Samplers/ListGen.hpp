/*
 * ListGen.hpp
 *
 *  Created on: Jan 23, 2018
 *      Author: guther
 */

#ifndef SRC_HEADERS_LISTGEN_HPP_
#define SRC_HEADERS_LISTGEN_HPP_

#include "Sampler.hpp"

namespace networkVMC{
template <typename F=std::complex<double>, typename coeffType=std::complex<double>>
class ListGen : public Sampler<F, coeffType>{
  public:
	ListGen(ExcitationGenerator const &eG_, Basis const &fullBasis_, detType const &HF,
			Parametrization<F, coeffType> const &para_, int numDets_=100);
	ListGen(Hamiltonian const &H_, Basis const &fullBasis_, detType const &HF,
			Parametrization<F, coeffType> const &para_, int numDets_=100);
	virtual ~ListGen();
	// create a dynamic polymorphic copy
	virtual ListGen<F, coeffType>* clone() const {return new ListGen<F, coeffType>(*this);}
	// get the i-th entry
	virtual void iterate(coeffType &cI, detType &dI, double& weight, int i);
	void diffuse(std::vector<detType> &list) const;
	void setDiffuseList(std::vector<detType > const &list){diffuseList=list;};
	virtual detType getDet(int i) const;
	virtual detType getDet() const;
    virtual int getNumDets() const;

    void resetSpecs();
private:
	mutable std::vector<detType > diffuseList;
  // sampling depends on the coefficients, as they have to be given alongside the determinants
    Parametrization<F, coeffType> const *para;
	mutable size_t pos;
    // and the corresponding basis including the information on the number of electrons
	// with a given spin
	Basis const *fullBasis;
    using Sampler<F, coeffType>::numDets;
};

}

#endif /* SRC_HEADERS_LISTGEN_HPP_ */
