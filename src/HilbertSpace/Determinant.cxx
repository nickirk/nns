/*
 * Determinant.cxx
 * Copyright (c) 2017
 * Author: Ke Liao and Kai Guther
 * All rights reserved
 */

#include "Determinant.hpp"

#include <vector>
#include <stdio.h>
#include <iostream>
#include <math.h>

namespace networkVMC{

void annihilate(detType &det, int pos){
  if (pos > static_cast<int>(det.size()) || pos < 0){
    std::cerr << "Error! Annihilate outside of range! "; 
    throw OutOfRangeError(pos);
  }
  if (det[pos]){ 
    det[pos] = 0;
  }
  else{ 
    std::cerr << "Error! Cannot annihilate an unoccupied state! ";
    throw InvalidAnnihilation(pos);
  }
}

//---------------------------------------------------------------------------------------------------//

void create(detType &det, int pos){
  if (pos > static_cast<int>(det.size()) || pos < 0){
    std::cerr << "Error! Create outside of range! "; 
    throw OutOfRangeError(pos);
  }
  if (det[pos]){
    std::cerr << "Error! Cannot create on an occupied state! ";
    throw InvalidCreation(pos);
  }
  else{
    det[pos] = 1;
  }
}

//---------------------------------------------------------------------------------------------------//

std::vector<int> getOccupiedPositions(detType const &det) {
  std::vector<int> positions;
  for (size_t i=0; i < det.size(); i++){
    if (det[i]) positions.push_back(i);
  }
  return positions;
}

//---------------------------------------------------------------------------------------------------//

int verbatimCast(detType const & det){
	int tmp = 0;
	for(size_t i = 0;i<det.size();++i){
		if(det[i]){
			tmp += pow(2,i);
		}
	}
	return tmp;
}

//---------------------------------------------------------------------------------------------------//

bool operator==(detType const& a, detType const& b) { return verbatimCast(a) == verbatimCast(b);}

bool operator < (detType const& lhs, detType const& rhs)
{
  return verbatimCast(lhs) < verbatimCast(rhs);
}

bool compare_det(detType const& lhs, detType const& rhs){

    return (lhs<rhs);
}

}