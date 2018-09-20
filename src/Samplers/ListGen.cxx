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
#include "../utilities/RNGWrapper.hpp"

namespace networkVMC{

ListGen::ListGen(ExcitationGenerator const &eG_, Basis const &fullBasis_, detType const &HF,
		Parametrization const &para_, int numDets_):
	Sampler(eG_,numDets_),para(&para_),pos(0),fullBasis(&fullBasis_){
	std::vector<detType> tmp(numDets_,HF);
	diffuse(tmp);
}

//---------------------------------------------------------------------------------------------------//

ListGen::ListGen(Hamiltonian const &H_, Basis const &fullBasis_, detType const &HF,
		Parametrization const &para_, int numDets_):
	Sampler(H_,HF,numDets_),para(&para_),pos(0),fullBasis(&fullBasis_){
	std::vector<detType> tmp(numDets_,HF);
	diffuse(tmp);
}

//---------------------------------------------------------------------------------------------------//

ListGen::~ListGen() {
}

//---------------------------------------------------------------------------------------------------//

void ListGen::iterate(coeffType &cI, detType &dI, double& weight, int i){
	// Fetch the next entry from the pre-arranged list
	dI = getDet(i);
	// Get its coefficient
	cI = para->getCoeff(dI);
}

//---------------------------------------------------------------------------------------------------//

detType ListGen::getDet(int i) const{
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
	  throw errors::OutOfRangeError(i);
  }
}

//---------------------------------------------------------------------------------------------------//

int ListGen::getNumDets() const{return diffuseList.size();}

//---------------------------------------------------------------------------------------------------//

void ListGen::diffuse(std::vector<detType> &list) const{
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
 RNGWrapper rng;

 for (size_t i = 0; i<list.size(); ++i){
  prandom = rng();
   c_i=para->getCoeff(list[i]);
   buf = Sampler::getRandomConnection(list[i],pEx);
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

}
