/*
 * FermiHubbardHamiltonian.cxx
 *
 *  Created on: Jun 20, 2018
 *      Author: guther
 */

#include "FermiHubbardHamiltonian.hpp"
#include "../HilbertSpace/Determinant.hpp"
#include "../utilities/Errors.hpp"

namespace networkVMC {

// the hubbard
double FermiHubbardHamiltonian::operator()(detType const &a, detType const &b) const{
	// first, see if a, b are valid
	sizeCheck(a.size(),2*(grid->size()));
	sizeCheck(a.size(),b.size());

	// then, get the excitation
	// holes/excits of b with respect to a
	std::vector<int> same, holes, excitations;
	// get the spins that differ between a and b
	getExcitation(a,b,excitations,holes,same);

	// RS hubbard is nearest-neighbours
	if(excitations.size() > 1) return 0.0;
	if(excitations.size() == 1){
		// the site indices of the hole and the site the electron moved to
		// why the @#$!?& are these 1-based
		// also, convert spin-site indices to spatial site-indices (lattices are always spatial)
		int hole = spatialIndex(holes[0] - 1);
		int part = spatialIndex(excitations[0] - 1);
		// only adjacent sites are coupled
		if(grid->isAdjacent(hole,part) ){
			// get the sign induces by the hopping
			auto sign = excitationSign(a,holes[0]-1,excitations[0]-1);
			// then return hopping matrix element
			return t*sign;
		}
		// all others have 0 matrix element
		return 0.0;
	}

	// if we reach this, a==b (or, equivalently, excitations.size() == 0)
	// That means, sum up all on-site contributions
	double diagElement = 0;
	// loop over all sites
	for(int i = 0; i<grid->size(); ++i){
		// check if a site is doubly occupied
		if(a[2*i] and a[2*i+1]){
			// if yes, add U to the diagonal
			diagElement += U;
		}
	}
	return diagElement;
}

//---------------------------------------------------------------------------------------------------//

void FermiHubbardHamiltonian::addCoupledStates(std::vector<detType> &list, detType const &source, int siteA, int siteB) const{
	// check for both spins if a hopping is possible
	for(int i = 0; i < 2; ++i){
		int indA = 2*siteA+i;
		int indB = 2*siteB+i;
		this->LatticeHamiltonian::addCoupledStates(list,source,indA,indB);
	}
}


} /* namespace networkVMC */
