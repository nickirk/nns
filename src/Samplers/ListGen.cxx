/*
 * ListGen.cxx
 *
 *  Created on: Jan 23, 2018
 *      Author: guther
 */

#include "ListGen.hpp"
#include "../HilbertSpace/Basis.hpp"
#include "../Network/Parametrization.hpp"
#include "../utilities/Errors.hpp"

namespace networkVMC{

template <typename T>
ListGen<T>::ListGen(ExcitationGenerator const &eG_, Basis const &fullBasis_, detType const &HF,
		Parametrization<T> const &para_, int numDets_):
	Sampler(eG_,HF,numDets_),para(&para_),pos(0),fullBasis(&fullBasis_){
	std::vector<detType> tmp(numDets_,HF);
	diffuse(tmp);
}

//---------------------------------------------------------------------------------------------------//

template <typename T>
ListGen<T>::ListGen(Hamiltonian const &H_, Basis const &fullBasis_, detType const &HF,
		Parametrization<T> const &para_, int numDets_):
	Sampler(H_,HF,numDets_),para(&para_),pos(0),fullBasis(&fullBasis_){
	std::vector<detType> tmp(numDets_,HF);
	diffuse(tmp);
}

//---------------------------------------------------------------------------------------------------//


template <typename T>
ListGen<T>::~ListGen() {
}

//---------------------------------------------------------------------------------------------------//

template <typename T>
void ListGen<T>::iterate(coeffType &cI, detType &dI, double& weight, int i){
	// Fetch the next entry from the pre-arranged list
	dI = getDet(i);
	// Get its coefficient
	cI = para->getCoeff(dI);
}

//---------------------------------------------------------------------------------------------------//

template <typename T>
detType ListGen<T>::getDet(int i) const{
  if(i < diffuseList.size() && i >= 0){
    return diffuseList[i];
  }
  else if(i >= diffuseList.size()){
	  // if we exceed the list size, generate a new one and start from the beginning
	  diffuse(diffuseList);
	  return getDet(0);
  }
  else{
	  // we passed a negative index, this is trouble
	  throw OutOfRangeError(i);
  }
}

//---------------------------------------------------------------------------------------------------//

template <typename T>
detType ListGen<T>::getDet() const{
	pos += 1;
	// cycle through the list, if we reach the end, start from the beginning
	if(pos > diffuseList.size()){
		pos = 1;
		// and also generate a new list
		diffuse(diffuseList);
	}
	return diffuseList[pos-1];
}

//---------------------------------------------------------------------------------------------------//

template <typename T>
int ListGen<T>::getNumDets() const{return diffuseList.size();}

//---------------------------------------------------------------------------------------------------//

template <typename T>
void ListGen<T>::diffuse(std::vector<detType> &list) const{
 list=diffuseList;
 detType buf;
 while (list.size() < static_cast<unsigned int>(numDets)){
   buf = getRandomDeterminant(*fullBasis);
   list.push_back(buf);
 }
 coeffType c_i = coeffType();
 coeffType c_j = coeffType();
 double prandom = 0.0;
 double pEx = 0.0;
 std::random_device rd;     // only used once to initialise (seed) engine
 std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
 std::uniform_real_distribution<double> uni;		// Uniform distribution from 0.0 to 1.0

 for (size_t i = 0; i<list.size(); ++i){
  prandom = uni(rng);
   c_i=para->getCoeff(list[i]);
   buf = getRandomConnection(list[i],pEx);
   //buf = getRandomDeterminant(spinConfig);
   c_j=para->getCoeff(buf);
   //getRandomCoupledState(buf,probUnbias);
   if (prandom - std::norm(c_j/c_i)<-1e-8){
     list[i]=buf;
   }
 }
 removeDuplicate(list);
 diffuseList = list;
}

//---------------------------------------------------------------------------------------------------//

template<typename T>
void ListGen<T>::resetSpecs(){
	// roll back the state to the initial one
	pos = 0;
}

//---------------------------------------------------------------------------//
//instantiate class
template class ListGen<VecType>;
template class ListGen<VecCType>;
}
