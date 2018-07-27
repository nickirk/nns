/*
 * HeisenbergHamiltonian.cxx
 *
 *  Created on: Jul 26, 2018
 *      Author: guther
 */

#include "HeisenbergHamiltonian.hpp"
#include "../HilbertSpace/Determinant.hpp"

namespace networkVMC {

double HeisenbergHamiltonian::operator ()(detType const &a, detType const &b) const{
	// the usual size-check
	if(a.size() != numSites or b.size() != numSites or b.size() != a.size()) throw SizeMismatchError(a.size(),numSites);
	// holes/excits of b with respect to a
	std::vector<int> same, holes, excitations;
	// get the spins that differ between a and b
	getExcitation(a,b,excitations,holes,same);
	// Heisenberg Hamiltonian only couples nearest neighbors
	if(excitations.size() > 1) return 0.0;
	// two possibilities: a==b or single excit
	if(excitations.size() == 1){
		// most likely: single excit
		// then: local excitations have -J as matrix element
		if(isAdjacent(excitations[0]-1,holes[0]-1)){
			// for some stupid reason, the excitations/holes
			// were set up to be 1-based, ugh
			return -J;
		}
		// all other have 0
		return 0.0;
	}

	// diagonal contribution: sum over all pairs of adjacent sites
	double diagElement = 0.0;
	for(int i = 0; i < numSites; ++i){
		// for each site, loop over the adjacent sites
		for(auto j:adjacencyList[i]){
			diagElement += msVal(a[i])*msVal(a[j]);
		}
	}
	// we counted each connection twice now, so use a factor J/2
	diagElement *= J/2.0;
	return diagElement;
}

//---------------------------------------------------------------------------------------------------//

bool HeisenbergHamiltonian::isAdjacent(int i, int j) const{
	// use the lookup table to determine if i, j are adjacent
	if(adjacencyList[i].find(j) != adjacencyList[i].end()) return true;
	return false;
}

//---------------------------------------------------------------------------------------------------//

// can be overloaded for different geometries
std::set<int> HeisenbergHamiltonian::adjacentSites(int i) const{
	// return all sites adjacent to i
	// use a set because each site can only appear once
	// lattice indices are 0-based

	// we map the sites canonically to integers
	// envision the n-d lattice as an n-dimensional tensor, then the sites are the entries
	// we label them as they would be addressed when stored as a contiguous array
	std::set<int> sites;
	int offset = 1;
	int tmp = 0;
	// loop over all dimensions and get the two additional sites adjacent in that dimension
	for(auto it = latticeDim.begin(); it != latticeDim.end(); ++it){
		// we can go right/left in the current dimension -> loop over direction
		for(int pre = -1; pre < 2; pre+=2){
			tmp = i + pre*offset;
			// check if we would go over the edge
			if(( i + (1+pre)/2*offset )%( offset*(*it) ) == 0 ){
				tmp = tmp - pre*offset*(*it);
			}
			// another way of going over the edge
			if(tmp < 0 or tmp >= numSites ){
				tmp = tmp - pre*offset*(*it);
			}
			// this is now the periodically continued adjacent site
			sites.insert(tmp);
		}
		offset *= *it;
	}

	return sites;
}

//---------------------------------------------------------------------------------------------------//

std::vector<detType> HeisenbergHamiltonian::getCoupledStates(detType const &source) const{
	// for each site, check if we can couple to an adjacent site
	std::vector<detType> AllCoupledStates;
	detType coupledState;
	for(int i = 0; i < numSites; ++i){
		for(auto j:adjacencyList[i]){
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


} /* namespace networkVMC */
