/*
 * LatticeHamiltonian.cxx
 *
 *  Created on: Jul 27, 2018
 *      Author: guther
 */

#include "../HilbertSpace/Determinant.hpp"
#include "LatticeHamiltonian.hpp"

namespace networkVMC{

std::vector<detType> LatticeHamiltonian::getCoupledStates(detType const &source) const{
	// for each site, check if we can couple to an adjacent site
	std::vector<detType> AllCoupledStates;
	detType coupledState;
	for(int i = 0; i < grid->size(); ++i){
		for(auto j:grid->adjacents(i)){
			// we can create a coupled det if and only if i,j have different spins
			// i<j to prevent double counting
			if( i< j and (source[i] xor source[j])){
				coupledState = source;
				// flip the spins at i,j
				if(source[i]){
					annihilate(coupledState,i);
					create(coupledState,j);
				}
				else{
					annihilate(coupledState,j);
					create(coupledState,i);
				}
				// we got a couple state, add it to the list
				AllCoupledStates.push_back(coupledState);
			}
		}
	}

	return AllCoupledStates;
}

}



