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
  detType bufPrev;
  int numRef=cDet.size();
  int random_integer(0);
  //std::cout << "size of list= " << list.size() << std::endl;
  //if(numStates%cDet.size()!=0) numComp = (numStates/cDet.size()+1)*cDet.size();
  //std::cout << "size of numComp= " << numComp << std::endl;
  if (numRef < numStates) {
    for (int i(0); i<numRef; ++i){
      list[i] = cDet[i];
    }
    std::random_device rng;     // only used once to initialise (seed) engine
    double const normalization=static_cast<double>(rng.max());
    random_integer=static_cast<int>(rng()/normalization*(numRef-1));
    buf = cDet[random_integer]; 
    bufPrev = buf;
    for (int i(numRef); i < numStates; ++i){
      while (verbatimCast(bufPrev)==verbatimCast(buf)){
        buf = getRandomDeterminant(buf);
      }
      list[i] = buf;
      bufPrev = buf;
    }
  }
  else if (numRef == numStates) list = cDet;
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


detType Sampler::getRandomDeterminant(detType const &startingPoint) const{
	double p{0};
	return getRandomCoupledState(startingPoint, p);
}

void removeDuplicate(std::vector<detType> &list){
 std::sort( list.begin(), list.end() );
 auto it=std::unique( list.begin(), list.end() );
 list.erase( it, list.end() );
}
