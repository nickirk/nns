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
    int random_integer;
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

detType Sampler::getRandomDeterminant(detType const &startingPoint) const{
  std::vector<detType > tempDets(0);
  // first, convert the determinant to an index
  int j = fullBasis.getIndexByDet(startingPoint);
  // get the range that corresponds to this index' row
  // here, we need the sorted representation. This is best done in the Hamiltonian class
  int lower = H.lowerPos(j);
  int upper = H.upperPos(j);
  std::random_device rd;
  std::random_device rng;
  int random_integer;
  std::mt19937 rng1(rd());    // random-number engine used (Mersenne-Twister in this case)
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
  int N=upper-lower; //number of connected determinants besides itself
  std::uniform_int_distribution<int> uni(lower, upper-1); // guaranteed unbiased
  random_integer = uni(rng1); 
  //std::cout << "N=" << N << std::endl;
  //std::cout << "random=" << random_integer << std::endl;
  while(tempDets.size() == 0){
      random_integer = uni(rng1); 
    //for(int i=lower;i<upper;++i){
      p=rng()/normalizer;
      // get all coupled determinants and accept them into the temporary list
      // with probability K_ij/K
      H.sparseAccess(random_integer,row,col,value);
      if (row == j) continue;
      //std::cout<<"Choosing "<<dblAbs(value)/K<<std::endl;
      //if(p<dblAbs(value)/K){
      if (p<1){
      //if(p<1./N){
	// here, we need the reverse of intcast, i.e. the conversion of the index
	// to the determinant. It shall just return lookuptable(i) for a given index i
	// (lookuptable contains the determinants for conversion to int)
	tempDets.push_back(fullBasis.getDetByIndex(row));
	//std::cout << "pushed back " << fullBasis.getIndexByDet(tempDets.back()) << std::endl;
      }
   //}
  }
  // pick a random determinant from the temporary list
  int const chosen=static_cast<int>(rng()/normalizer*tempDets.size());
  return tempDets[chosen];
}

void Sampler::removeDuplicate(std::vector<detType> &list){
 std::sort( list.begin(), list.end() );
 list.erase( std::unique( list.begin(), list.end() ), list.end() );
}
