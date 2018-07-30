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
#include <vector>

namespace networkVMC {

LatticeExcitgen::~LatticeExcitgen() {
}

//---------------------------------------------------------------------------------------------------//

detType LatticeExcitgen::generateExcitation(detType const& source,
		double& pGen) {
	// random number generator
	RNGWrapper rng();
	// first, pick a random site in source
	// that is, get the occupied sites
	auto occSites = getOccupiedPositions(source);
	// and pick one at random
	int sourceSite = occSites[rng()*occSites.size()];
	// accumulate 1/pGen and then do pGen = 1.0/pGen
	pGen = occSites.size();
	// now, pick a random adjacent site that is empty
	auto selects = getOccAdjs(sourceSite, source, HL);
	// if we picked invalid, return the source
	if(selects.size()==0) return source;
	// assign the correct pGen
	pGen = 1.0/(pGen*selects.size());
	// now pick a random coupled site
	int targetSite = selects[rng()*selects.size()];
	detType target = source;
	// apply the excitation
	annihilate(target,sourceSite);
	// if all went well, return the excited det
	return target;
}

//---------------------------------------------------------------------------------------------------//

double LatticeExcitgen::getExcitationProb(detType const& source,
		detType const& target) {
	auto occSites = getOccupiedPositions(source);
	// get the exctiation
	std::vector<int> holes, excits, same;
	getExcitation(source,target,excits,holes,same);
	// probability of picking the hole
	double pGen = occSites.size();
	// consistency check, move to debug later
	if(holes.size() != excits.size()) throw SizeMismatchError(holes.size(),excits.size());
	// check if there is an excitation at all
	if(holes.size() != 0){
		// probability of picking the excit
		auto selects = getOccAdjs(holes[0], source, HL);
		// if there is a hole, there has to be an excitation
		pGen *= selects.size();
	}

}

//---------------------------------------------------------------------------------------------------//

std::vector<int> getOccAdjs(int i, detType const &source, LatticeHamiltonian const &HL){
	auto adjs = HL.adjacents(i);
	std::vector<int> selects(0);
	for(auto site:adjs){
		// only empty sites count
		if(!source[site]){
			selects.push_back(site);
		}
	}
	return selects;
}

} /* namespace networkVMC */
