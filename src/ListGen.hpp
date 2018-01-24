/*
 * ListGen.hpp
 *
 *  Created on: Jan 23, 2018
 *      Author: guther
 */

#ifndef SRC_LISTGEN_HPP_
#define SRC_LISTGEN_HPP_

#include "Nnw.hpp"
#include "Determinant.hpp"
#include "Sampler.hpp"

class ListGen : public Sampler{
public:
	ListGen(Hamiltonian const &H_, Basis const &fullBasis_, int numDets_, detType const &HF, NeuralNetwork const &NNW);
	virtual ~ListGen();
	void iterate(coeffType &cI, detType &dI, NeuralNetwork const &NNW) const;
	void diffuse(std::vector<detType> &list, std::vector<int> const& spinConfig);
	detType getDet() const;
private:
	NeuralNetwork const &NNW;
	std::vector<detType > diffuseList;
	mutable size_t pos;
};

#endif /* SRC_LISTGEN_HPP_ */
