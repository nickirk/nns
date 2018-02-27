/*
 * Sampler.cxx
 * Copyright (c) 2017
 * Author: Ke Liao and Kai Guther
 * All rights reserved
 */
#include <vector>
#include <random>
#include <algorithm>
#include <Eigen/Dense>
#include "Sampler.hpp"
#include "State.hpp"
#include "Determinant.hpp"
#include <iostream>

void Sampler::initialiseList(std::vector<detType> &list, std::vector<int> const& spinConfig){
  list.clear();
  detType buf;
  std::random_device rngd;
  int numRef(cDet.size());
  if (numRef < numDets) {
    for (int i(0); i<numRef; ++i){
      list.push_back(cDet[i]);
    }
  //else{ throw an error: ref cannot be larger than numDets
  }
    //std::random_device rng;     // only used once to initialise (seed) engine
    //double const normalization=static_cast<double>(rng.max());
    //random_integer=static_cast<int>(rng()/normalization*(numRef-1));
    //for (int i(numRef); i < numDets; ++i){
  while (list.size() < numDets){
   buf = getRandomDeterminant(spinConfig);
   list.push_back(buf);
  }
}

void Sampler::diffuse(std::vector<detType> &list, std::vector<int> const& spinConfig){
 removeDuplicate(list);
 detType buf;
 //std::cout << "before adding random dets=" << list.size() << std::endl;
 while (list.size() < numDets){
   buf = getRandomDeterminant(spinConfig);
   list.push_back(buf);
 }
 //std::cout << "after adding random dets=" << list.size() << std::endl;
 coeffType c_i;
 coeffType c_j;
 Eigen::VectorXd lastLayerActivation;
 std::random_device rngd;
 double const normalizerd = static_cast<double>(rngd.max());
 double prandom=rngd()/normalizerd;
 double probGet(0.);
 double probUnbias(0.);
 for (size_t i(0); i<list.size(); ++i){
   prandom=rngd()/normalizerd;
   lastLayerActivation=NNW.feedForward(list[i]);
   c_i=coeffType(lastLayerActivation[0], lastLayerActivation[1]);
   buf = H.getRandomCoupledState(list[i],probUnbias);
   //buf = getRandomDeterminant(spinConfig);
   lastLayerActivation = NNW.feedForward(buf);
   c_j=coeffType(lastLayerActivation[0], lastLayerActivation[1]);
   //getRandomCoupledState(buf,probUnbias);
   if (prandom - std::pow(std::norm(c_j),2)/std::pow(std::norm(c_i),2)<-1e-8){
     list[i]=buf;
   }
 }
 //std::cout << "after diffuse=" << list.size() << std::endl;
 removeDuplicate(list);
 //std::cout << "after remove dups=" << list.size() << std::endl;
}
/*
void Sampler::generateList(std::vector<detType > &list) const{
  // this just repeatedly gets new random states and adds them to the list
  // if the number of target states is not an integer multiple of the number of
  // reference states, we have to round up the number of target states
  //int numComp = numDets;
  //list = std::vector<detType >(numDets);
  list.clear();
  detType buf;
  detType bufPrev;
  int numRef=cDet.size();
  int random_integer(0);
  //if(numDets%cDet.size()!=0) numComp = (numDets/cDet.size()+1)*cDet.size();
  //std::cout << "size of numComp= " << numComp << std::endl;
  std::random_device rd;     // only used once to initialise (seed) engine
  std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
  std::uniform_int_distribution<int> uni(0,numRef-1); // guaranteed unbiased
  std::random_device rngd;
  double const normalizerd = static_cast<double>(rngd.max());
  if (numRef < numDets) {
    for (int i(0); i<numRef; ++i){
      list.push_back(cDet[i]);
    }
    //std::random_device rng;     // only used once to initialise (seed) engine
    //double const normalization=static_cast<double>(rng.max());
    //random_integer=static_cast<int>(rng()/normalization*(numRef-1));
    //for (int i(numRef); i < numDets; ++i){
   double prandom=rngd()/normalizerd;
   if (prandom-0.8 < -1e-8){
    random_integer = uni(rng);
    buf = cDet[random_integer]; 
   }
   else{
    buf = getRandomDeterminant(3,3,16);
   }
    bufPrev = buf;
    while (list.size() < numDets){
      prandom=rngd()/normalizerd;
      if (prandom-1. < -1e-8){
      //if (true){
        //for (int depth(0);  depth < 4; depth++){ 
          while (verbatimCast(bufPrev)==verbatimCast(buf) ){
            buf = bufPrev;
            buf = getRandomConnection(buf);
          }
          list.push_back(buf);
          bufPrev = buf;
        //}
      }
      else {
        buf = getRandomDeterminant(3,3,16);
        list.push_back(buf);
      }
    }
  }
  else if (numRef == numDets) list = cDet;
}
*/

detType Sampler::getRandomConnection(detType const &startingPoint) const{
	double p{0};
	return H.getRandomCoupledState(startingPoint, p);
}

void removeDuplicate(std::vector<detType> &list){
 std::sort( list.begin(), list.end() );
 auto it=std::unique( list.begin(), list.end() );
 list.erase( it, list.end() );
}
/*
void getFullSpace(std::vector<int> const &list, std::vector<int> const &spinConfig){
  int numSpinUp(spinConfig[0]);
  int numSpinDown(spinConfig[1]);
  int numStates(spinConfig[2]);
  detType buf(numStates, 0);  
  //int numSpinUp = (numEle%2) ? numEle/2 : (numEle+1)/2;
  std::vector<int> possibleSpinUpStates;
  for (int i(0);i<(numStates/2); ++i){
    possibleSpinUpStates.push_back(i*2);
  }
  //int numSpinDown = numEle - numSpinUp;
  std::vector<int> possibleSpinDownStates;
  for (int i(0);i<(numStates/2); ++i){
    possibleSpinDownStates.push_back(i*2+1);
  }
  for (int i(0); i<numSpinUp; ++i){
    create(randomDet,possibleSpinUpStates[random_integer]);
    possibleSpinUpStates.erase(possibleSpinUpStates.begin()+random_integer);
  }
  for (int i(0); i<numSpinDown; ++i){
    std::uniform_int_distribution<int> uni(0,possibleSpinDownStates.size()-1); // guaranteed unbiased
    random_integer=uni(rng);
    create(randomDet,possibleSpinDownStates[random_integer]);
    possibleSpinDownStates.erase(possibleSpinDownStates.begin()+random_integer);
  }
}
*/
detType getRandomDeterminant(std::vector<int> const &spinConfig){ 
  int numSpinUp(spinConfig[0]);
  int numSpinDown(spinConfig[1]);
  int numStates(spinConfig[2]);
  detType randomDet(numStates, 0);  
  //int numSpinUp = (numEle%2) ? numEle/2 : (numEle+1)/2;
  std::vector<int> possibleSpinUpStates;
  for (int i(0);i<(numStates/2); ++i){
    possibleSpinUpStates.push_back(i*2);
  }
  //int numSpinDown = numEle - numSpinUp;
  std::vector<int> possibleSpinDownStates;
  for (int i(0);i<(numStates/2); ++i){
    possibleSpinDownStates.push_back(i*2+1);
  }
  int random_integer(0);
  //std::vector<int> selectedSpinUpStates;
  std::random_device rd;     // only used once to initialise (seed) engine
  std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
  std::uniform_int_distribution<int> uni(0,possibleSpinUpStates.size()-1); // guaranteed unbiased
  for (int i(0); i<numSpinUp; ++i){
    //std::random_device rng;     // only used once to initialise (seed) engine
    //double const normalization=static_cast<double>(rng.max());
    //random_integer=static_cast<int>(rng()/normalization*(possibleSpinUpStates.size()-1));
    //std::cout << "size =" << possibleSpinUpStates.size() << std::endl;
    //selectedSpinUpStates.push_back(possibleSpinUpStates[random_integer]); 
    std::uniform_int_distribution<int> uni(0,possibleSpinUpStates.size()-1); // guaranteed unbiased
    random_integer=uni(rng);
    create(randomDet,possibleSpinUpStates[random_integer]);
    //std::cout << "random int =" << random_integer << std::endl;
    //std::cout << "possibleSpinUpStates[random_integer] =" << possibleSpinUpStates[random_integer]<< std::endl;
    //std::cout << "intCast after create=" <<verbatimCast(randomDet) <<std::endl;
    possibleSpinUpStates.erase(possibleSpinUpStates.begin()+random_integer);
  }
  //std::vector<int> selectedSpinDownStates;
  for (int i(0); i<numSpinDown; ++i){
    std::uniform_int_distribution<int> uni(0,possibleSpinDownStates.size()-1); // guaranteed unbiased
    random_integer=uni(rng);
    //selectedSpinDownStates.push_back(possibleSpinDownStates[random_integer]); 
    create(randomDet,possibleSpinDownStates[random_integer]);
    //std::cout << "size =" << possibleSpinDownStates.size() << std::endl;
    //std::cout << "random int =" << random_integer << std::endl;
    //std::cout << "possibleSpinDownStates[random_integer] =" << possibleSpinDownStates[random_integer]<< std::endl;
    //std::cout << "intCast after create=" <<verbatimCast(randomDet) <<std::endl;
    possibleSpinDownStates.erase(possibleSpinDownStates.begin()+random_integer);
  }
  return randomDet; 
}
