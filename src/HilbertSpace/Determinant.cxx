/*
 * Determinant.cxx
 * Copyright (c) 2017
 * Author: Ke Liao and Kai Guther
 * All rights reserved
 */

#include "Determinant.hpp"
#include "../utilities/Errors.hpp"
#include <vector>
#include <stdio.h>
#include <iostream>
#include <math.h>

namespace networkVMC{

void annihilate(detType &det, int pos){
  // remove an electron from orbital pos from det
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
  // add an electron to orbital pos from det
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
  // return the occupied orbital's indices
  std::vector<int> positions;
  for (size_t i=0; i < det.size(); i++){
    if (det[i]) positions.push_back(i);
  }
  return positions;
}

//---------------------------------------------------------------------------------------------------//

int verbatimCast(detType const & det){
	// 1:1 conversion of a determinant to integer (not surjective)
	int tmp = 0;
	for(size_t i = 0;i<det.size();++i){
		if(det[i]){
			tmp += pow(2,i);
		}
	}
	return tmp;
}

//---------------------------------------------------------------------------------------------------//

int getExcitLvl(detType const &a, detType const &b){
	// excit lvls are only defined between same-sized dets
	if(a.size() != b.size()) throw SizeMismatchError(a.size(),b.size());
	// with the same number of electrons
	if(getOccupiedPositions(a).size() != getOccupiedPositions(b).size())
		throw InvalidDeterminantError();
	int elvl{0};
	for(std::size_t i=0; i < a.size(); ++i){
		// if the dets differ, this counts towards the excitation lvl
		if(a[i]^b[i]) elvl += 1;
	}
	return elvl/2;
}

void getExcitation(detType const &a, detType const &b, std::vector<int> &excitations,
		std::vector<int> &holes, std::vector<int> &same){
	int diff = 0;
	int const d = a.size();
	// we can only get excitations between same-sized determinants
	if(d != b.size()) throw SizeMismatchError(d,b.size());
    // get the differences in occupied obitals
    for (int i=0;i<d;++i){
        diff=static_cast<int>(a[i])-static_cast<int>(b[i]);
        if (diff>0){
            // holes created in alpha
            holes.push_back(i+1);
        }
        if (diff<0){
            // particles created in beta
            excitations.push_back(i+1);
        }
        else{
            // these spin orbitals are present in both
            if (diff==0 && a[i]){
                same.push_back(i+1);
            }
        }
    }
}

//---------------------------------------------------------------------------------------------------//

// binary operators for determinants

bool operator==(detType const& a, detType const& b) { return verbatimCast(a) == verbatimCast(b);}

bool operator < (detType const& lhs, detType const& rhs)
{
  return verbatimCast(lhs) < verbatimCast(rhs);
}

bool compare_det(detType const& lhs, detType const& rhs){

    return (lhs<rhs);
}

// for debugging: output a determinant

void printDet(detType const &out){
	auto orbs = getOccupiedPositions(out);
	std::cout<<"(";
	for(int i : orbs){
		std::cout << i << ", ";
	}
	std::cout<<")"<<std::endl;
}

}
