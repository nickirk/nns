/*
 * ListGen.cxx
 *
 *  Created on: Jan 23, 2018
 *      Author: guther
 */

#include "ListGen.hpp"
#include "../Network/Parametrization.hpp"

namespace networkVMC{

ListGen::ListGen(Hamiltonian const &H_, Basis const &fullBasis_,
		detType const &HF, Parametrization const &para_, int numDets_):
	Sampler(H_,fullBasis_,HF,para_,numDets_),pos(0){
	std::vector<detType> tmp(numDets_,HF);
	diffuse(tmp);
}

//---------------------------------------------------------------------------------------------------//

ListGen::~ListGen() {
}

//---------------------------------------------------------------------------------------------------//

void ListGen::iterate(coeffType &cI, detType &dI) const{
	// Fetch the next entry from the pre-arranged list
	dI = getDet();
	// Get its coefficient
	cI = para->getCoeff(dI);
}

//---------------------------------------------------------------------------------------------------//

detType ListGen::getDet(int i) const{
  return diffuseList[i];
}

//---------------------------------------------------------------------------------------------------//

detType ListGen::getDet() const{
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

 std::random_device rd;     // only used once to initialise (seed) engine
 std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
 std::uniform_real_distribution<double> uni;		// Uniform distribution from 0.0 to 1.0

 for (size_t i = 0; i<list.size(); ++i){
  prandom = uni(rng);
   c_i=para->getCoeff(list[i]);
   buf = getRandomConnection(list[i]);
   //buf = getRandomDeterminant(spinConfig);
   c_j=para->getCoeff(buf);
   //getRandomCoupledState(buf,probUnbias);
   if (prandom - std::pow(std::norm(c_j),2)/std::pow(std::norm(c_i),2)<-1e-8){
     list[i]=buf;
   }
 }
 removeDuplicate(list);
 diffuseList = list;
}

}
