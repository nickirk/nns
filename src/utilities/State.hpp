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

// A state is a list of determinants and the corresponding list of
// their coefficients and weights, also the lists of the respective coupled Dets and
// coupled coeffs. We want the State to behave like a std::vector.
class State{
public:
	State():
    storedDets(std::vector<detType>(0)), 
    storedCoeffs(std::vector<coeffType>(0)),
    storedWeights(std::vector<double>(0)), 
    fSortedCoeff(0), fSortedDet(0), fSortedWeight(0){};
	State(int size_):State(){resize(size_);}
	State(detType const &det_, coeffType const &coeff_):
	storedDets(std::vector<detType>(1,det_)),
    storedCoeffs(std::vector<coeffType>(1,coeff_)),
    storedWeights(std::vector<double>(1,1)),
    fSortedCoeff(0), fSortedDet(0), fSortedWeight(0){};
	State( std::vector<detType> const &dets_, 
         std::vector<coeffType> const &coeffs_
         ):
    storedDets(dets_),storedCoeffs(coeffs_), storedWeights(coeffs_.size(),1),
    fSortedCoeff(0), fSortedDet(0), fSortedWeight(0){
		// A state has to have one coefficient per determinant
		// one might argue that supplying less coefficient should be fine and that
		// the rest should be filled with zeroes, we might change that
	    if ( storedDets.size() != storedCoeffs.size() && 
           storedCoeffs.size() != storedWeights.size() &&
           storedWeights.size() != storedDets.size() ) 
        throw SizeMismatchError(storedDets.size(),storedCoeffs.size());
	};

//---------------------------------------------------------------------------------------------------//

	// access operators for determinants and coefficients of the state
	detType const& det(int i) const{return storedDets[i];}
	coeffType const& coeff(int i) const{return storedCoeffs[i];}
	double const& weight(int i) const{return storedWeights[i];}

	// implement the non-const version via the const version
	detType& det(int i) {
    return const_cast<detType &>(static_cast<State const &>((*this)).det(i));
  }
	coeffType& coeff(int i){
    return const_cast<coeffType &>(static_cast<State const &>((*this)).coeff(i));
  }
	double& weight(int i){
    return const_cast<double &>(static_cast<State const &>((*this)).weight(i));
  }

	// access operators for the coupled determinants
	std::vector<detType>const& coupledDets(int i) const{return cDets[i];}
	std::vector<coeffType>const& coupledCoeffs(int i) const{return cCoeffs[i];}
	std::vector<double>const& coupledWeights(int i) const{return cWeights[i];}

	// same as above: prevent code duplication
	std::vector<detType>& coupledDets(int i) {
    return const_cast<std::vector<detType > &>(
        static_cast<State const &>((*this)).coupledDets(i)
        );
  }
	std::vector<coeffType>& coupledCoeffs(int i){
    return const_cast<std::vector<coeffType> &>(
        static_cast<State const &>((*this)).coupledCoeffs(i)
        );
  }
	std::vector<double>& coupledWeights(int i){
    return const_cast<std::vector<double> &>(
        static_cast<State const &>((*this)).coupledWeights(i)
        );
  }

	// vector utilities
	void clear() {
    storedDets.clear(); storedCoeffs.clear(); storedWeights.clear(); 
    cDets.clear(); cCoeffs.clear(); cWeights.clear();
  }
	void resize(int i) {
    storedDets.resize(i); storedCoeffs.resize(i); storedWeights.resize(i);
    cDets.resize(i); cCoeffs.resize(i); cWeights.resize(i);
  }

	// vector utilities
	size_t size()const{return storedDets.size();}

//---------------------------------------------------------------------------------------------------//


//---------------------------------------------------------------------------------------------------//
	// sort the storedDets + storedCoeffs according to determinant order
	void sortDet(){
	  // first, get the permutation used for sorting
	  auto perm = getPermutation(storedDets,detComparer());
	  // then sort the dets and coeffs according to it (and cDets/cCoeffs with them)
	  permuteAll(perm);
    fSortedDet = 1;
    fSortedCoeff = 0;
    fSortedWeight = 0;
	}
//---------------------------------------------------------------------------------------------------//
	// or according to highest coefficient
  void sortCoeffs(){
  // sort the state using the coefficient-wise comparison
  // we employ the sorting.hpp utilities (same as above)
    auto perm = getPermutation(storedCoeffs,coeffComparer());
    permuteAll(perm);
    fSortedCoeff = 1;
    fSortedDet = 0;
    fSortedWeight = 0;
  }

//---------------------------------------------------------------------------------------------------//
  void sortWeights(){
    auto perm = getPermutation(storedWeights, weightComparer());
    permuteAll(perm);
    fSortedWeight = 1;
    fSortedDet = 0;
    fSortedCoeff = 0;
  }

//---------------------------------------------------------------------------------------------------//
    // this is just to prevent code duplication
    void permuteAll(std::vector<std::size_t> const &perm){
    // reorder the content of the state according to the permutation perm
      storedDets = applyPermutation(storedDets,perm);
      storedCoeffs =  applyPermutation(storedCoeffs,perm);
      cDets = applyPermutation(cDets,perm);
      cCoeffs = applyPermutation(cCoeffs,perm);
      storedWeights = applyPermutation(storedWeights, perm);
    }
 //totalSize is numDets+num of all coupled Dets
 int totalSize() const{
   int size(0);
   for (size_t i(0); i<storedDets.size();++i){
     size+=1;
     for (size_t j(0); j < cDets[i].size(); ++j){
       size+=1;
     }
   }
   return size;
 };

 size_t locate(size_t iDet) const{
  //locate the i'th det in memory
  // TODO: Directly use cDets[i].size() - the current loop makes
  // no sense: We call cDets[i].size() to compute cDets[i].size()???
   size_t size(0);
   for (size_t i(0); i< iDet;++i){
     for (size_t j(0); j < cDets[i].size(); ++j){
       size+=1;
     }
   }
   return size;
 };
//---------------------------------------------------------------------------------------------------//
private:
  // whenever a new det is added, update its coeff and weight
	std::vector<detType> storedDets;
	std::vector<coeffType> storedCoeffs;
    std::vector<double> storedWeights;
   
  // util of flags
    bool fSortedCoeff;
    bool fSortedDet;
    bool fSortedWeight;


	// preliminary implementation of coupled stuff
  // whenever a new det is added, update its coeff and weight
	std::vector<std::vector<detType> > cDets;
	std::vector<std::vector<coeffType >> cCoeffs;
    std::vector<std::vector<double>> cWeights;
};

}

#endif /* SRC_UTILITIES_STATE_HPP_ */
