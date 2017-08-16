/*
 * Sampler.cxx
 * Copyright (c) 2017
 * Author: Ke Liao and Kai Guther
 * All rights reserved
 */
#include <vector>
#include <random>
#include "Hamiltonian.hpp"
#include "Sampler.hpp"

void sampler::generateList(std::vector<detType > &list) const{
  // this just repeatedly gets new random states and adds them to the list
  // if the number of target states is not an integer multiple of the number of
  // reference states, we have to round up the number of target states
  if(numStates%cDet.size()!=0) numStates = (numStates/cDet.size()+1)*cDet.size()
  list = std::vector<detType >(numStates);
  detType buf;
  for(size_t j=0;j<cDet.size();++j){
    // now, numStates/cDet.size() is always an integer
    buf = cDet[j];
    for(int i=0;i<numStates/cDet.size();++i){
      buf=getRandomDeterminant(buf);
      list[i]=buf;
    }
  }
}

detType sampler::getRandomDeterminant(detType const &startingPoint) const{
  std::vector<detType > tempDets(0);
  // first, convert the determinant to an index
  int j = fullBasis.getIndexByDet(startingPoint);
  // get the range that corresponds to this index' row
  // here, we need the sorted representation. This is best done in the Hamiltonian class
  int lower = H.lowerPos(j);
  int upper = H.upperPos(j);
  int row,col;
  std::random_device rng;
  double const normalizer = static_cast<double>(rng.max());
  double value,p;
  double K=0.0;
  // compute the normalization sum_j |K_ij|
  for(int i=lower;i<=upper;++i){
    H.sparseAccess(i,row,col,value);
    K+=dblAbs(value);
  }
  for(int i=lower;i<=upper;++i){
    p=rng()/normalizer;
    // get all coupled determinants and accept them into the temporary list
    // with probability K_ij/K
    H.sparseAccess(i,row,col,value);
    if(p<value/K){
      // here, we need the reverse of intcast, i.e. the conversion of the index
      // to the determinant. It shall just return lookuptable(i) for a given index i
      // (lookuptable contains the determinants for conversion to int)
      tempDets.push_back(fullBasis.getDetByIndex(col));
    }
  }
  // pick a random determinant from the temporary list
  int const chosen=static_cast<int>(rng()/normalizer*tempDets.size());
  return tempDets[chosen];
}
