/*
 * ListGen.hpp
 *
 *  Created on: Jan 23, 2018
 *      Author: guther
 */

#ifndef SRC_HEADERS_LISTGEN_HPP_
#define SRC_HEADERS_LISTGEN_HPP_

#include "../HilbertSpace/Determinant.hpp"
#include "../Network/Nnw.hpp"
#include "../utilities/SpinConfig.hpp"
#include "Sampler.hpp"

namespace networkVMC{

class ListGen : public Sampler{
public:
	ListGen(Hamiltonian const &H_, Basis const &fullBasis_, int numDets_, detType const &HF, NeuralNetwork const &NNW);
	virtual ~ListGen();
	virtual void iterate(coeffType &cI, detType &dI) const;
	void diffuse(std::vector<detType> &list) const;
	void setDiffuseList(std::vector<detType > const &list){diffuseList=list;};
	virtual detType getDet(int i) const;
	virtual detType getDet() const;
    virtual int getNumDets() const;
private:
	NeuralNetwork const &NNW;
	mutable std::vector<detType > diffuseList;
	mutable size_t pos;
};

}

#endif /* SRC_HEADERS_LISTGEN_HPP_ */
