/*
 * HeisenbergHamiltonian.cxx
 *
 *  Created on: Jul 26, 2018
 *      Author: guther
 */

#include "Lattice.hpp"
#include "../HilbertSpace/Determinant.hpp"

namespace networkVMC {

//---------------------------------------------------------------------------------------------------//

bool Lattice::isAdjacent(int i, int j) const{
	// use the lookup table to determine if i, j are adjacent
	if(adjacencyList[i].find(j) != adjacencyList[i].end()) return true;
	return false;
}

//---------------------------------------------------------------------------------------------------//

// can be overloaded for different geometries
std::set<int> Lattice::adjacentSites(int i) const{
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
			if(pbc or tmp == (i + pre*offset)) sites.insert(tmp);
		}
		offset *= *it;
	}

	return sites;
}


} /* namespace networkVMC */
