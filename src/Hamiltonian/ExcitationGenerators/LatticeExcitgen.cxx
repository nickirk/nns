/*
 * LatticeExcitgen.cxx
 *
 *  Created on: Jul 30, 2018
 *      Author: guther
 */

#include "LatticeExcitgen.hpp"
#include "../../utilities/RNGWrapper.hpp"
#include "../../utilities/Errors.hpp"
#include "../../HilbertSpace/Determinant.hpp"
#include "../LocalQN.hpp"
#include <vector>
#include <iostream>

namespace networkVMC {

LatticeExcitgen::~LatticeExcitgen() {
}

//---------------------------------------------------------------------------------------------------//

detType LatticeExcitgen::generateExcitation(detType const& source,
		double& pGen) {
	// random number generator
	RNGWrapper rng;
	// first, pick a random site in source
	// that is, get the occupied sites
	auto occSites = getOccupiedPositions(source);
	// and pick one at random
	int sourceSite = occSites[rng()*occSites.size()];
	// accumulate 1/pGen and then do pGen = 1.0/pGen
	pGen = occSites.size();
	// now, pick a random adjacent site that is empty
	// with the correct local QN
	LocalQN lQN(sourceSite,mapTypeToRange(HL.type()));
	// adjacency is spatial
	auto selects = getOccAdjs(sourceSite, source, HL,lQN);
	// if we picked invalid, return the source
	if(selects.size()==0){
		// and make sure to set the right pGen
		pGen = 1.0/pGen;
		return source;
	}
	// assign the correct pGen
	pGen = 1.0/(pGen*selects.size());
	// now pick a random coupled site
	int targetSite = selects[rng()*selects.size()];
	detType target = source;
	// apply the excitation
	annihilate(target, sourceSite);
	create(target, targetSite);
	// if all went well, return the excited det
	return target;
}

//---------------------------------------------------------------------------------------------------//

double LatticeExcitgen::getExcitationProb(detType const& source,
		detType const& target) {
	auto occSites = getOccupiedPositions(source);
	// probability of picking the hole
	double pGen = occSites.size();
	// get the exctiation
	std::vector<int> holes, excits, same;
	getExcitation(source,target,excits,holes,same);
	// consistency check, move to debug later
	if(holes.size() != excits.size()) throw errors::SizeMismatchError(holes.size(),excits.size());

	// check if there is an excitation at all
	if(holes.size() != 0){
		LocalQN lQN(holes[0]-1,mapTypeToRange(HL.type()));
		// probability of picking the excit
		auto selects = getOccAdjs(holes[0]-1, source, HL, lQN);
		// if there is a hole, there has to be an excitation
		if(selects.size()==0){
			std::cout << "ERROR, unable to get excit from " << holes[0]-1 <<std::endl;
			printDet(source);
			printDet(target);
		}
		pGen *= selects.size();
	}
	return 1.0/pGen;
}

//---------------------------------------------------------------------------------------------------//

std::vector<int> getOccAdjs(int i, detType const &source, LatticeHamiltonian const &HL, LocalQN const &lQN){
	auto adjs = HL.adjacents(lQN.spatialIndex(i));
	std::vector<int> selects(0);
	for(auto site:adjs){
		// only empty sites count
		// now, we need the orbital indices
		if(!source[lQN.orbitalIndex(site)]){
			selects.push_back(lQN.orbitalIndex(site));
		}
	}
	return selects;
}

} /* namespace networkVMC */
