/*
 * ListGen.cxx
 *
 *  Created on: Jan 23, 2018
 *      Author: guther
 */

#include "Sampler.hpp"
#include "ListGen.hpp"

ListGen::ListGen(Hamiltonian const &H_, Basis const &fullBasis_, int numDets_, detType const &HF, NeuralNetwork const &NNW_):
	Sampler(H_,fullBasis_,HF,numDets_),NNW(NNW_),pos(0){
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
	cI = NNW.getCoeff(dI);
	// Dont forget to cache the state
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
 while (list.size() < numDets){
   buf = getRandomDeterminant(sC);
   list.push_back(buf);
 }
 coeffType c_i;
 coeffType c_j;
 Eigen::VectorXd lastLayerActivation;
 std::random_device rngd;
 double const normalizerd = static_cast<double>(rngd.max());
 double prandom=rngd()/normalizerd;
 double probUnbias(0.);
 for (size_t i(0); i<list.size(); ++i){
   prandom=rngd()/normalizerd;
   c_i=NNW.getCoeff(list[i]);
   buf = getRandomConnection(list[i]);
   //buf = getRandomDeterminant(spinConfig);
   c_j=NNW.getCoeff(buf);
   //getRandomCoupledState(buf,probUnbias);
   if (prandom - std::pow(std::norm(c_j),2)/std::pow(std::norm(c_i),2)<-1e-8){
     list[i]=buf;
   }
 }
 removeDuplicate(list);
 diffuseList = list;
}
