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
template<typename T=VecType>
class ListGen : public Sampler{
public:
	ListGen(Hamiltonian const &H_, Basis const &fullBasis_, detType const &HF,
			Parametrization<T> const &para_, int numDets_=100);
	virtual ~ListGen();
	virtual void iterate(coeffType &cI, detType &dI) const;
	void diffuse(std::vector<detType> &list) const;
	void setDiffuseList(std::vector<detType > const &list){diffuseList=list;};
	virtual detType getDet(int i) const;
	virtual detType getDet() const;
    virtual int getNumDets() const;
private:
	mutable std::vector<detType > diffuseList;
  // sampling depends on the coefficients, as they have to be given alongside the determinants
    Parametrization<T> const *para;
	mutable size_t pos;
};

}

#endif /* SRC_HEADERS_LISTGEN_HPP_ */
