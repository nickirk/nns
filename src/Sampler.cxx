/*
 * Sampler.cxx
 * Copyright (c) 2017
 * Author: Ke Liao and Kai Guther
 * All rights reserved
 */
#include <vector>
#include <random>
#include "Sampler.hpp"
#include <iostream>

void Sampler::generateList(std::vector<detType > &list) const{
  // this just repeatedly gets new random states and adds them to the list
  // if the number of target states is not an integer multiple of the number of
  // reference states, we have to round up the number of target states
  int numComp = numStates;
  std::cout << "size of list= " << list.size() << std::endl;
  if(numStates%cDet.size()!=0) numComp = (numStates/cDet.size()+1)*cDet.size();
  std::cout << "size of numComp= " << numComp << std::endl;
  list = std::vector<detType >(numComp);
  detType buf;
  int numRef=cDet.size();
  std::cout << "size of cDet= " << numRef << std::endl;
  int numBatch=numComp/numRef;
  for(size_t j=0;j<numRef;++j){
    // now, numStates/cDet.size() is always an integer
    buf = cDet[j];
    list[j] = buf;
    for(int i=0;i<numBatch-1;++i){
      buf=getRandomDeterminant(buf);
      list[numRef+j*(numBatch-1)+i]=buf;
    }
  }
}

detType Sampler::getRandomDeterminant(detType const &startingPoint) const{
  std::vector<detType > tempDets(0);
  // first, convert the determinant to an index
  int j = fullBasis.getIndexByDet(startingPoint);
  // get the range that corresponds to this index' row
  // here, we need the sorted representation. This is best done in the Hamiltonian class
  int lower = H.lowerPos(j);
  int upper = H.upperPos(j);
  std::random_device rng;
  int row = -1;
  int col = -1;
  double const normalizer = static_cast<double>(rng.max());
  double value,p;
  double K=0.0;
  // compute the normalization sum_j |K_ij|
  for(int i=lower;i<upper;++i){
    H.sparseAccess(i,row,col,value);
    if (row == j) continue;
    K+=dblAbs(value);
  }
  while(tempDets.size() == 0){
    for(int i=lower;i<upper;++i){
      p=rng()/normalizer;
      // get all coupled determinants and accept them into the temporary list
      // with probability K_ij/K
      H.sparseAccess(i,row,col,value);
      if (row == j) continue;
      //std::cout<<"Choosing "<<dblAbs(value)/K<<std::endl;
      if(p<dblAbs(value)/K){
	// here, we need the reverse of intcast, i.e. the conversion of the index
	// to the determinant. It shall just return lookuptable(i) for a given index i
	// (lookuptable contains the determinants for conversion to int)
	tempDets.push_back(fullBasis.getDetByIndex(row));
	//std::cout << "pushed back " << fullBasis.getIndexByDet(tempDets.back()) << std::endl;
      }
    }
  }
  // pick a random determinant from the temporary list
  int const chosen=static_cast<int>(rng()/normalizer*tempDets.size());
  return tempDets[chosen];
}
