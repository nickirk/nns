/*
 * Sampler.cxx
 * Copyright (c) 2017
 * Author: Ke Liao and Kai Guther
 * All rights reserved
 */
#include <vector>
#include <random>
#include <algorithm>
#include "Sampler.hpp"
#include "Determinant.hpp"
#include <iostream>

void Sampler::generateList(std::vector<detType > &list) const{
  // this just repeatedly gets new random states and adds them to the list
  // if the number of target states is not an integer multiple of the number of
  // reference states, we have to round up the number of target states
  //int numComp = numStates;
  list = std::vector<detType >(numStates);
  detType buf;
  int numRef=cDet.size();
  std::random_device rd;     // only used once to initialise (seed) engine
  //std::cout << "size of list= " << list.size() << std::endl;
  //if(numStates%cDet.size()!=0) numComp = (numStates/cDet.size()+1)*cDet.size();
  //std::cout << "size of numComp= " << numComp << std::endl;
  if (numRef < numStates) {
    for (int i(0); i<numRef; ++i){
      list[i] = cDet[i];
    }
    int random_integer{0};
    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
    std::uniform_int_distribution<int> uni(0,numRef-1); // guaranteed unbiased
    random_integer = uni(rng); 
    buf = cDet[random_integer]; 
    for (int i(numRef); i < numStates; ++i){
      //std::cout << "random_integer= " << random_integer <<std::endl;
      //std::cout << "size of Ref= " << numRef <<std::endl;
      //std::cout << "size of Dets= " << numStates <<std::endl;
      //std::cout << "i= " << i <<std::endl;
      //std::cout << "intcast ref= " << verbatimCast(buf) << std::endl;
      buf = getRandomDeterminant(buf);
      list[i] = buf;
      //std::cout << "intcast new= " << verbatimCast(buf) << std::endl;
    }
  }
  else if (numRef == numStates) list = cDet;
  //std::cout << "size of cDet= " << numRef << std::endl;
  //int numBatch=numComp/numRef;
  //for(size_t j=0;j<numRef;++j){
  //  // now, numStates/cDet.size() is always an integer
  //  buf = cDet[j];
  //  list[j] = buf;
  //  for(int i=0;i<numBatch-1;++i){
  //    buf=getRandomDeterminant(buf);
  //    list[numRef+j*(numBatch-1)+i]=buf;
  //  }
  //}
}

void Sampler::removeDuplicate(std::vector<detType> &list){
 std::sort( list.begin(), list.end() );
 list.erase( std::unique( list.begin(), list.end() ), list.end() );
}

detType Sampler::getRandomDeterminant(detType const &startingPoint) const{
	double p{0};
	return getRandomCoupledState(startingPoint, p);
}
