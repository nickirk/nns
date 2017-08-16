/*
 * Determinant.cxx
 * Copyright (c) 2017
 * Author: Ke Liao and Kai Guther
 * All rights reserved
 */

#include <vector>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include "Determinant.hpp"

void annihilate(detType &det, int pos){
  if (pos > det.size() || pos < 0){
    std::cerr << "Error! Annihilate outside of range! "; 
    throw outOfRangeError(pos);
  }
  if (det[pos]){ 
    det[pos] = 0;
  }
  else{ 
    std::cerr << "Error! Cannot annihilate an unoccupied state! ";
    throw invalidAnnihilation(pos);
  }
}

void create(detType &det, int pos){
  if (pos > det.size() || pos < 0){
    std::cerr << "Error! Create outside of range! "; 
    throw outOfRangeError(pos);
  }
  if (det[pos]){
    std::cerr << "Error! Cannot create on an occupied state! ";
    throw invalidCreation(pos);
  }
  else{
    det[pos] = 1;
  }
}

std::vector<int> getOccupiedPositions(detType const &det) {
  std::vector<int> positions;
  for (int i=0; i < det.size(); i++){
    if (det[i]) positions.push_back(i);
  }
  return positions;
}

int verbatimCast(detType const & det){
	int tmp = 0;
	for(size_t i = 0;i<det.size();++i){
		if(det[i]){
			tmp += pow(2,i);
		}
	}
	return tmp;
}

bool operator==(detType const& a, detType const& b) { return verbatimCast(a) == verbatimCast(b);}
