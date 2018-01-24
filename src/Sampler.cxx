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

detType Sampler::getRandomConnection(detType const &startingPoint) const{
	double p{0};
	return getRandomCoupledState(startingPoint, p);
}

void removeDuplicate(std::vector<detType> &list){
 std::sort( list.begin(), list.end() );
 auto it=std::unique( list.begin(), list.end() );
 list.erase( it, list.end() );
}

detType getRandomDeterminant(std::vector<int> const &spinConfig){ 
  int numSpinUp(spinConfig[0]);
  int numSpinDown(spinConfig[1]);
  int numStates(spinConfig[2]);
  detType randomDet(numStates, 0);  
  std::vector<int> possibleSpinUpStates;
  for (int i(0);i<(numStates/2); ++i){
    possibleSpinUpStates.push_back(i*2);
  }
  std::vector<int> possibleSpinDownStates;
  for (int i(0);i<(numStates/2); ++i){
    possibleSpinDownStates.push_back(i*2+1);
  }
  int random_integer(0);
  std::random_device rd;     // only used once to initialise (seed) engine
  std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
  std::uniform_int_distribution<int> uni(0,possibleSpinUpStates.size()-1); // guaranteed unbiased
  for (int i(0); i<numSpinUp; ++i){

    //selectedSpinUpStates.push_back(possibleSpinUpStates[random_integer]); 
    std::uniform_int_distribution<int> uni(0,possibleSpinUpStates.size()-1); // guaranteed unbiased
    random_integer=uni(rng);
    create(randomDet,possibleSpinUpStates[random_integer]);

    possibleSpinUpStates.erase(possibleSpinUpStates.begin()+random_integer);
  }
  for (int i(0); i<numSpinDown; ++i){
    std::uniform_int_distribution<int> uni(0,possibleSpinDownStates.size()-1); // guaranteed unbiased
    random_integer=uni(rng);
    //selectedSpinDownStates.push_back(possibleSpinDownStates[random_integer]); 
    create(randomDet,possibleSpinDownStates[random_integer]);
    possibleSpinDownStates.erase(possibleSpinDownStates.begin()+random_integer);
  }
  return randomDet; 
}
