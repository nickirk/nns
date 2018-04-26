/*
 * State.hpp
 *
 *  Created on: Oct 25, 2017
 *      Author: guther
 */

#ifndef SRC_UTILITIES_STATE_HPP_
#define SRC_UTILITIES_STATE_HPP_

#include <vector>

#include "../HilbertSpace/Determinant.hpp"
#include "TypeDefine.hpp"
#include "Errors.hpp"

namespace networkVMC{

// A state is a list of determinants and their coefficients
class State{
public:
	State():dets(std::vector<detType>(0)), coeffs(std::vector<coeffType>(0)){};
	State(detType const &det_, coeffType const &coeff_):
		dets(std::vector<detType>(1,det_)), coeffs(std::vector<coeffType>(1,coeff_)) {};
	State(std::vector<detType> const &dets_, std::vector<coeffType> const &coeffs_):dets(dets_),coeffs(coeffs_){
		// A state has to have one coefficient per determinant
		// one might argue that supplying less coefficient should be fine and that
		// the rest should be filled with zeroes, we might change that
	    if(dets.size() != coeffs.size()) throw SizeMismatchError(dets.size(),coeffs.size());
	};

	// access operators for determinants and coefficients of the state
	detType const& getDet(int i) const{return dets[i];}
	coeffType const& getCoeff(int i) const{return coeffs[i];}

	// implement the non-const version via the const version
	detType& getDet(int i) {return const_cast<detType &>(static_cast<State const &>((*this)).getDet(i));}
	coeffType& getCoeff(int i){return const_cast<coeffType &>(static_cast<State const &>((*this)).getCoeff(i));}

	// access operators for the coupled determinants
	std::vector<detType>const& coupledDets(int i) const{return cdets[i];}
	std::vector<coeffType>const& coupledCoeffs(int i) const{return ccoeffs[i];}

	// same as above: prevent code duplication
	std::vector<detType>& coupledDets(int i) {return const_cast<std::vector<detType > &>(static_cast<State const &>((*this)).coupledDets(i));}
	std::vector<coeffType>& coupledCoeffs(int i){return const_cast<std::vector<coeffType> &>(static_cast<State const &>((*this)).coupledCoeffs(i));}

	// vector utilities
	void clear() {dets.clear(); coeffs.clear(); cdets.clear(); ccoeffs.clear();}
	void resize(int i) {dets.resize(i); coeffs.resize(i); cdets.resize(i); ccoeffs.resize(i);}
	size_t size()const{return dets.size();}
private:
	std::vector<detType> dets;
	std::vector<coeffType> coeffs;

	// preliminary implementation of coupled stuff
	std::vector<std::vector<detType> > cdets;
	std::vector<std::vector<coeffType >> ccoeffs;
};

}

#endif /* SRC_UTILITIES_STATE_HPP_ */
