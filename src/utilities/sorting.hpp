/*
 * sorting.hpp
 *
 *  Created on: Apr 27, 2018
 *      Author: guther
 */

#ifndef SRC_UTILITIES_SORTING_HPP_
#define SRC_UTILITIES_SORTING_HPP_

#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include "TypeDefine.hpp"
#include "../HilbertSpace/Determinant.hpp"

namespace networkVMC{

// get the permutation (as a vector of ints) that is used to sort sortVec
// with a functor comp()
template <typename T, typename compare>
std::vector<std::size_t> getPermutation(std::vector<T> const &sortVec,
		compare const &comp){
// initially, there is no permutation, so set the permutation to a vector
// of numerically increasing ints
  std::vector<std::size_t> permutation(sortVec.size());
  std::iota(permutation.begin(),permutation.end(),0);
// here comes the trick: we sort permutation, but as a comparison,
// we use the comparison between elements of sortVec
  std::sort(permutation.begin(),permutation.end(),
	[&](std::size_t i, std::size_t j){return comp(sortVec[i],sortVec[j]);});
	// for this, inline the comparison as lambda
  return permutation;
}

//---------------------------------------------------------------------------------------------------//

// reorder a given vector sortVec according to some permutation
template <typename T>
std::vector<T> applyPermutation(
		std::vector<T> &sortVec,
		std::vector<std::size_t> const &permutation
		){
// the output shall be sortVec
  std::vector<T> newVec(sortVec.size());
// but reordered according to permutation
// for this, use std::transform: it applies a given function to the range
// from permutation.begin() to permutation.end() and writes into newVec
  std::transform(permutation.begin(),permutation.end(),newVec.begin(),
    [&](std::size_t i){return sortVec[i];});
  // the function we give just returns the entries of sortVec
  // at the position given by the permutation
  return newVec;
}

//---------------------------------------------------------------------------------------------------//

// functor for comparing dets
class detComparer{
  public:
  detComparer(){};
  // this is a callable object which does nothing but return if a<b for dets
  bool operator()(detType const &a, detType const &b) const{return a>b;}
};

//---------------------------------------------------------------------------------------------------//

// functor for comparing coeffs
template <typename coeffType=std::complex<double>>
class coeffComparer{
  public:
  coeffComparer(){};
  // this is another callable object which only returns if a<b (in the sense of absolute value)
  bool operator()(coeffType const &a, coeffType const &b) const{return std::abs(a)>std::abs(b);}
};

// functor for comparing weights
class weightComparer{
  public:
  weightComparer(){};
  // this is another callable object which only returns if a>b (in the sense of value)
  bool operator()(double const &a, double const &b) const{return a>b;}
};

}

#endif /* SRC_UTILITIES_SORTING_HPP_ */
