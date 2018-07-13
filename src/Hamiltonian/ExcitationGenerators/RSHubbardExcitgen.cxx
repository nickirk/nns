/*
 * RSHubbardExcitgen.cxx
 *
 *  Created on: Jun 19, 2018
 *      Author: guther
 */

#include "RSHubbardExcitgen.hpp"

namespace networkVMC {

RSHubbardExcitgen::RSHubbardExcitgen():ClonableExcitgen<RSHubbardExcitgen>() {
}

RSHubbardExcitgen::~RSHubbardExcitgen() {
}

detType RSHubbardExcitgen::generateExcitation(
		detType const &source,  double &pGet){
  std::random_device rng;
  double const normalization=static_cast<double>(rng.max());
  int const d=source.size();
  detType target=source;
  //for hubbard
  //hubbard one-particle basis: 1(up), 1(down), 2(up), ..., L(up), L(down)
  std::vector<int> spawnLeft, spawnRight;
  // generate a list of all possible spawns
  getRSHubSpawnLists(source, spawnLeft, spawnRight);

  int const numSpawn=spawnLeft.size()+spawnRight.size();
  // We cannot do any excitations from this determinant
  if(numSpawn == 0){
	  pGet = 1.0;
	  return source;
  }
  int p = static_cast<int>(rng()/normalization*numSpawn);
  //pick one of those sites at random (including direction)
  pGet=1.0/static_cast<double>(numSpawn);
  int const offset=spawnLeft.size();
  int newPos = 0;
  if(p>=offset){
    newPos=spawnRight[p-offset];
    annihilate(target,newPos);
    create(target,(newPos+2)%d);
  }
  else{
    newPos=spawnLeft[p];
    annihilate(target,newPos);
    create(target,(newPos-2+d)%d);
  }
  return target;
}

//---------------------------------------------------------------------------------------------------//

double RSHubbardExcitgen::getExcitationProb(
			detType const &source, detType const &target){
	// the excitation probability is 1/(number of possible hoppings)
	std::vector<int> spawnLeft, spawnRight;
	// so first, get a list of all possible hoppings
	getRSHubSpawnLists(source, spawnLeft, spawnRight);

	// and return the inverse
	int numSpawns = spawnLeft.size() + spawnRight.size();
	return 1.0/static_cast<double>(numSpawns);
}

//---------------------------------------------------------------------------------------------------//

void getRSHubSpawnLists(detType const &source, std::vector<int> &spawnLeft,
		std::vector<int> &spawnRight){
	  int d = source.size();
	  for(int i=0;i<d;++i){
	    //search all sites from which one can hop to the left/right
	    if(source[i]){
	      if(source[(i+2)%d]==0){
		spawnRight.push_back(i);
	      }
	      if(source[(i-2+d)%d]==0){
		spawnLeft.push_back(i);
	      }
	    }
	  }
}

} /* namespace networkVMC */
