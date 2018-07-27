/*
 * LatticeHamiltonian.cxx
 *
 *  Created on: Jul 27, 2018
 *      Author: guther
 */

#include "../HilbertSpace/Determinant.hpp"
#include "LatticeHamiltonian.hpp"

namespace networkVMC{

// this implementation works for lattices with one-to-one correspondence between single-particle states and sites
std::vector<detType> LatticeHamiltonian::getCoupledStates(detType const &source) const{
	// for each site, check if we can couple to an adjacent site
	std::vector<detType> AllCoupledStates;
	detType coupledState;
	for(int i = 0; i < grid->size(); ++i){
		for(auto j:grid->adjacents(i)){
			// we can create a coupled det if and only if i,j have different spins
			// i<j to prevent double counting
			if( i< j){
				addCoupledStates(AllCoupledStates,source,i,j);
			}
		}
	}
	return AllCoupledStates;
}

//---------------------------------------------------------------------------------------------------//

void LatticeHamiltonian::addCoupledStates(std::vector<detType> &list, detType const &source, int siteA, int siteB) const{
	// siteA and siteB are two adjacent sites: add all states reachable from source
	// by hopping on those sites
	// if we can hop, add the excitation
	if(source[siteA] xor source[siteB]){
		list.push_back(excite(source,siteA,siteB));
	}
}

}



