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
#include <algorithm>
#include <iostream>
#include <math.h>

namespace networkVMC{

void annihilate(detType &det, int pos){
  // remove an electron from orbital pos from det
  if (pos > static_cast<int>(det.size()) || pos < 0){
    std::cerr << "Error! Annihilate outside of range! "; 
    throw errors::OutOfRangeError(pos);
  }
  if (det[pos]){ 
    det[pos] = 0;
  }
  else{ 
    std::cerr << "Error! Cannot annihilate an unoccupied state! ";
    throw errors::InvalidAnnihilation(pos);
  }
}

//---------------------------------------------------------------------------------------------------//

void create(detType &det, int pos){
  // add an electron to orbital pos from det
  if (pos > static_cast<int>(det.size()) || pos < 0){
    std::cerr << "Error! Create outside of range! "; 
    throw errors::OutOfRangeError(pos);
  }
  if (det[pos]){
    std::cerr << "Error! Cannot create on an occupied state! ";
    throw errors::InvalidCreation(pos);
  }
  else{
    det[pos] = 1;
  }
}

//---------------------------------------------------------------------------------------------------//

int JWStringLength(detType const &a, int annihilatorIndex, int creatorIndex){
	  //construct sign for conversion canonical shape
	  //(basisState 1(up),1(down),...,L(up),L(down)) to relative shape (a_exc^\dagger a_holes source)
	  int fermiSign{0};
	  int start{0}, end{0};
	  fermiSign = 0;
	  if(annihilatorIndex>creatorIndex){
	    start=creatorIndex;
	    end=annihilatorIndex;
	  }
	  else{
	    start=annihilatorIndex;
	    end=creatorIndex;
	  }
	  for(int k=start+1;k<end;++k){
	    if(a[k]){
	      fermiSign += 1;
	    }
	  }
	  return fermiSign;
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
	if(a.size() != b.size()) throw errors::SizeMismatchError(a.size(),b.size());
	// with the same number of electrons
	if(getOccupiedPositions(a).size() != getOccupiedPositions(b).size())
		throw errors::InvalidDeterminantError(a);
	int elvl{0};
	for(std::size_t i=0; i < a.size(); ++i){
		// if the dets differ, this counts towards the excitation lvl
		if(a[i]^b[i]) elvl += 1;
	}
	return elvl/2;
}

//---------------------------------------------------------------------------------------------------//

void getExcitation(detType const &a, detType const &b, std::vector<int> &excitations,
		std::vector<int> &holes, std::vector<int> &same){
	int diff = 0;
	std::size_t const d = a.size();
	holes.clear();
	excitations.clear();
	same.clear();
	// we can only get excitations between same-sized determinants
	if(d != b.size()) throw errors::SizeMismatchError(d,b.size());
    // get the differences in occupied obitals
    for (std::size_t i=0;i<d;++i){
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

detType excite(detType const &source, int i, int j){
	// take source and excite i<->j (depending on which one is empty)
    detType coupledState = source;
	if(source[i]){
		annihilate(coupledState,i);
		create(coupledState,j);
	}
	else{
		annihilate(coupledState,j);
		create(coupledState,i);
	}
	return coupledState;
}

//---------------------------------------------------------------------------------------------------//

// binary operators for determinants

bool operator==(detType const& a, detType const& b) { return verbatimCast(a) == verbatimCast(b);}

bool operator< (detType const& lhs, detType const& rhs)
{
  return verbatimCast(lhs) < verbatimCast(rhs);
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

//---------------------------------------------------------------------------//

// Go through a list of dets and remove duplicates (generic function)
void removeDuplicate(std::vector<detType> &list){
 std::sort( list.begin(), list.end() );
// make use of std::unique for a sorted list
 auto it = std::unique( list.begin(), list.end() );
 list.erase( it, list.end() );
}

}
