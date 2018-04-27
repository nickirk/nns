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
#include "sorting.hpp"
#include "Errors.hpp"

namespace networkVMC{

// A state is a list of determinants and their coefficients
class State{
public:
	State():storedDets(std::vector<detType>(0)), storedCoeffs(std::vector<coeffType>(0)){};
	State(detType const &det_, coeffType const &coeff_):
		storedDets(std::vector<detType>(1,det_)), storedCoeffs(std::vector<coeffType>(1,coeff_)) {};
	State(std::vector<detType> const &dets_, std::vector<coeffType> const &coeffs_):storedDets(dets_),storedCoeffs(coeffs_){
		// A state has to have one coefficient per determinant
		// one might argue that supplying less coefficient should be fine and that
		// the rest should be filled with zeroes, we might change that
	    if(storedDets.size() != storedCoeffs.size()) throw SizeMismatchError(storedDets.size(),storedCoeffs.size());
	};

//---------------------------------------------------------------------------------------------------//

	// access operators for determinants and coefficients of the state
	detType const& det(int i) const{return storedDets[i];}
	coeffType const& coeff(int i) const{return storedCoeffs[i];}

	// implement the non-const version via the const version
	detType& det(int i) {return const_cast<detType &>(static_cast<State const &>((*this)).det(i));}
	coeffType& coeff(int i){return const_cast<coeffType &>(static_cast<State const &>((*this)).coeff(i));}

	// access operators for the coupled determinants
	std::vector<detType>const& coupledDets(int i) const{return cdets[i];}
	std::vector<coeffType>const& coupledCoeffs(int i) const{return ccoeffs[i];}

	// same as above: prevent code duplication
	std::vector<detType>& coupledDets(int i) {return const_cast<std::vector<detType > &>(static_cast<State const &>((*this)).coupledDets(i));}
	std::vector<coeffType>& coupledCoeffs(int i){return const_cast<std::vector<coeffType> &>(static_cast<State const &>((*this)).coupledCoeffs(i));}

	// vector utilities
	void clear() {storedDets.clear(); storedCoeffs.clear(); cdets.clear(); ccoeffs.clear();}
	void resize(int i) {storedDets.resize(i); storedCoeffs.resize(i); cdets.resize(i); ccoeffs.resize(i);}
	size_t size()const{return storedDets.size();}

//---------------------------------------------------------------------------------------------------//
	// sort the storedDets + storedCoeffs according to determinant order
	void sortDet(){
	// first, get the permutation used for sorting
	auto perm = getPermutation(storedDets,detComparer());
	// then sort the dets and coeffs according to it (and cdets/ccoeffs with them)
	permuteAll(perm);
	}
//---------------------------------------------------------------------------------------------------//
	// or according to highest coefficient
    void sortCoeffs(){
    // sort the state using the coefficient-wise comparison
    // we employ the sorting.hpp utilities (same as above)
      auto perm = getPermutation(storedCoeffs,coeffComparer());
      permuteAll(perm);
    }

//---------------------------------------------------------------------------------------------------//
    // this is just to prevent code duplication
    void permuteAll(std::vector<std::size_t> const &perm){
    // reorder the content of the state according to the permutation perm
      applyPermutation(storedDets,perm);
      applyPermutation(storedCoeffs,perm);
      applyPermutation(cdets,perm);
      applyPermutation(ccoeffs,perm);
    }

//---------------------------------------------------------------------------------------------------//
private:
	std::vector<detType> storedDets;
	std::vector<coeffType> storedCoeffs;

	// preliminary implementation of coupled stuff
	std::vector<std::vector<detType> > cdets;
	std::vector<std::vector<coeffType >> ccoeffs;
};

}

#endif /* SRC_UTILITIES_STATE_HPP_ */
