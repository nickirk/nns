/*
 * HeisenbergHamiltonian.cxx
 *
 *  Created on: Jul 27, 2018
 *      Author: guther
 */

#include "HeisenbergHamiltonian.hpp"
#include "../utilities/Errors.hpp"
#include "../HilbertSpace/Determinant.hpp"

namespace networkVMC {

double HeisenbergHamiltonian::operator ()(detType const &a, detType const &b) const{
	// the usual size-check
	sizeCheck(a.size(),grid->size());
	sizeCheck(b.size(),a.size());
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
		if(grid->isAdjacent(excitations[0]-1,holes[0]-1)){
			// for some stupid reason, the excitations/holes
			// were set up to be 1-based, ugh
			return -J;
		}
		// all other have 0
		return 0.0;
	}

	// diagonal contribution: sum over all pairs of adjacent sites
	double diagElement = 0.0;
	for(int i = 0; i < grid->size(); ++i){
		// for each site, loop over the adjacent sites
		for(auto j:grid->adjacents(i)){
			diagElement += msVal(a[i])*msVal(a[j]);
		}
	}
	// we counted each connection twice now, so use a factor J/2
	diagElement *= J/2.0;
	return diagElement;
}

} /* namespace networkVMC */
