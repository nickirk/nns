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

namespace networkVMC{

// A state contains a determinant and its coefficient, and might further contain information
// on coupled determinants
class State{
public:
	State():det(detType({0})), coeff(coeffType(0.,0.)), 
         coupledCoeffs(std::vector<coeffType>(0)), 
         coupledDets(std::vector<detType>(0)) {};
	State(detType const &det_, coeffType const &coeff_):
		det(det_), coeff(coeff_),
        coupledCoeffs(std::vector<coeffType>(0)),
        coupledDets(std::vector<detType>(0)) {};
	State(
              detType  const &det_, 
              coeffType const &coeff_, 
              std::vector<detType> const &coupledDets_,
              std::vector<coeffType> const &coupledCoeffs_
              ):
		det(det_), 
                coeff(coeff_), 
                coupledCoeffs(coupledCoeffs_),
                coupledDets(coupledDets_){
		// A state has to have one coefficient per determinant
		// one might argue that supplying less coefficient should be fine and that
		// the rest should be filled with zeroes, we might change that
	        // if(dets.size() != coeffs.size()) throw sizeMismatchError(dets.size(),coeffs.size());
	};
	//detType getDet(int i) const{return dets[i];}
	//coeffType getCoeff(int i) const{return coeffs[i];}
	//std::vector<coeffType> getAllCoeff() const{return coeffs;}
	//std::vector<detType> getDets() const {return dets;}
        //std::vector<coeffType> getCoupledCoeffs(int i) const{return coupledCoeffs[i];}
        //std::vector<detType> getCoupledDets(int i) const{return coupledDets[i];}
	//size_t size()const{return dets.size();}

	detType det;
	coeffType coeff;
// potentially also have a number of coupled determinants stored
        std::vector<coeffType> coupledCoeffs;
        std::vector<detType> coupledDets;
        bool operator < (State const &m) const {
	  return det < m.det;
	}
};

}

#endif /* SRC_UTILITIES_STATE_HPP_ */
