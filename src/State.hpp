/*
 * State.hpp
 *
 *  Created on: Oct 25, 2017
 *      Author: guther
 */

#ifndef SRC_STATE_HPP_
#define SRC_STATE_HPP_

#include <vector>
#include "Determinant.hpp"
#include "CoeffType.hpp"

class State{
public:
	State():dets(std::vector<detType >(0)), coeffs(std::vector<coeffType >(0)) {};
	State(std::vector<detType > const &dets_, std::vector<coeffType > const &coeffs_):
		dets(dets_), coeffs(coeffs_){
		// A state has to have one coefficient per determinant
		// one might argue that supplying less coefficient should be fine and that
		// the rest should be filled with zeroes, we might change that
		if(dets.size() != coeffs.size()) throw sizeMismatchError(dets.size(),coeffs.size());
	};
	detType getDet(int i) const{return dets[i];}
	coeffType getCoeff(int i) const{return coeffs[i];}
	size_t size()const{return dets.size();}
private:
	std::vector<detType > dets;
	std::vector<coeffType > coeffs;
};

#endif /* SRC_STATE_HPP_ */
